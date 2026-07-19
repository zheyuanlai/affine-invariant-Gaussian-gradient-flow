# Basepoint covariance and the projected-ray collision barrier

## 0. Verdict

Let \(\mu\) be isotropic and log-concave, and suppose that a median
signed-distance witness admits a normal-cell disintegration

\[
 \mu=\int \nu_q\,d\eta(q),\qquad
 X=Y_q+T_qN_q,\qquad |N_q|=1,                           \tag{0.1}
\]

where \(Y_q\) and \(N_q\) are quotient-measurable and \(0\) is a median
of every conditional law of \(T_q\).  The first two steps of the proposed
projection route are valid, with no dimension loss.

1. If \(m_q=\mathbb E_qT_q\) and
   \(\sigma_q^2=\operatorname {Var}_qT_q\), then

   \[
       |m_q|\le \sigma_q,\qquad
       \boxed{\operatorname {Cov}_\eta(Y_q)\preceq4I},
       \qquad
       \boxed{\mathbb E_\eta Y_qY_q^T\preceq4I}.       \tag{0.2}
   \]

   Log-concavity is not needed for these inequalities.  The proof is the
   law of total covariance plus the one-sided Chebyshev inequality.

2. Write

   \[
      D=\mathbb E|T_q|,
      \qquad K=D_1(\mu),
      \qquad D\ge\kappa K                              \tag{0.3}
   \]

   for a \(\kappa\)-near-maximal witness.  If \(C_M\) is the universal
   constant in E. Milman's implication

   \[
                         C_P(\mu)\le C_M D_1(\mu)^2,    \tag{0.4}
   \]

   then a set of quotient mass at least \(\kappa^2/(16C_M)\) satisfies

   \[
      {D\over2\sqrt2}\le\sigma_q
       \le {4C_M\over\kappa^2}D.                       \tag{0.5}
   \]

   Thus a fixed amount of the quotient really does live at the bad scale.

For balanced normal cells the transport-orientation tensor is

\[
 M=\int a_qN_qN_q^T\,d\eta(q),
 \qquad a_q=\mathbb E_q|T_q|,                           \tag{0.6}
\]

and it obeys

\[
       \operatorname {tr}M=D,
       \qquad \|M\|_{HS}\le1.                          \tag{0.7}
\]

Consequently its effective rank is **at least** \(D^2\).  This is the
opposite of a low-dimensional capture theorem.  From (0.7) one cannot
choose a rank-\(O(D^2)\) subspace which carries a fixed fraction of the
direction mass.  The algebraic obstruction is

\[
                     M={D\over N}I_N,qquad N\gg D^2.   \tag{0.8}
\]

It satisfies every assertion in (0.7), while every rank-\(cD^2\)
projection captures only \(cD^2/N=o(1)\) of its trace.  The isotropic
radial exponential law realizes this diffuse matrix exactly at its natural
constant scale.  Hence the proposed choice of \(H\) is not a consequence
of basepoint covariance, cell-scale localization, or the transport tensor.

There is nevertheless a sharp conditional result.  If an
\(O(D^2)\)-dimensional subspace \(H\) captures a fixed-mass packet of the
long-ray directions, then exactly one of the following occurs.

* Opposite projected tail packets remain separated by \(cD\).  They give a
  1-Lipschitz witness of size \(cD\) for the isotropic log-concave marginal
  \(P_H\mu\).  The known \(\sqrt{\log d}\) KLS bound in
  \(d=O(D^2)\) dimensions then forces \(D=O(1)\).
* A fixed amount of opposite tail mass can be coupled at projected distance
  at most \(cD\).  The corresponding full points are not close: the eikonal
  inequality forces their omitted-coordinate separation to be of order
  \(D\).  After removing the small omitted normal components, this produces
  a new transverse **basepoint** direction law of effective rank
  \(\Omega(D^2)\).  It produces no local perimeter or first-moment saving.

The second branch is the exact obstruction.  Projection has discarded the
large physical separation needed by any chord-splicing competitor.  The
covariance lemma merely permits the transverse displacement when it is
spread over \(\Omega(D^2)\) directions.  Iterating this observation can use
\(n/D^2\) mutually transverse packets, so it gives a depth depending on
\(n\), which is inadmissible for KLS.

The model audit is consistent and sharp.  A cube halfspace is the captured
affine branch.  A radial exponential median sphere is the diffuse-collision
branch and is semistable.  The atomic regular simplex saturates dimension
\(D^2\) and is a globally maximizing nonlinear shell coloring, but is not
log-concave.  Filling the simplex destroys cell balance and creates a fixed
contact contribution.  The product-exponential maximum is log-concave and
diffuse, but its regular cells are unbalanced and a fixed ridge bevel saves
perimeter.  These tests identify the missing assertion precisely: a true
large-scale maximizer would need an extremality-and-log-concavity theorem
which turns transverse collision rank into physical contact.  Neither
Paouris nor thin shell supplies that theorem.

## 1. Setting and hypotheses

All assertions are intrinsic to the affine support.  After whitening that
support, assume it is \(\mathbb R^n\) and

\[
                         \mathbb EX=0,qquad
                         \mathbb EXX^T=I.              \tag{1.1}
\]

Let \(Q\) be a standard Borel quotient and suppose (0.1) is a regular
conditional disintegration.  Thus, conditionally on \(q\), the only random
quantity in (0.1) is the scalar \(T_q\).  We assume

\[
 \mathbb P_q(T_q\ge0)\ge\frac12,
 \qquad
 \mathbb P_q(T_q\le0)\ge\frac12.                       \tag{1.2}
\]

For an absolutely continuous conditional law, exact normal-cell balance
means equality in both parts of (1.2).  Define

\[
 m_q=\mathbb E_qT_q,qquad
 \sigma_q^2=\mathbb E_q(T_q-m_q)^2,qquad
 B_q=\mathbb E[X\mid q]=Y_q+m_qN_q.                    \tag{1.3}
\]

Every quantity is square-integrable as a consequence of the identities
proved below.  No preliminary integrability of \(Y_q\) needs to be assumed.

For the scale theorem, assume additionally that the conditional laws are
one-dimensional log-concave and that

\[
                            f(X)=T_q                    \tag{1.4}
\]

for a globally 1-Lipschitz median-centered signed-distance witness.  A
normal Voronoi cell of a smooth stationary half-mass signed-distance
interface has exactly this form: its first variation bisects the cell, and
the normal coordinate equals signed distance until the cut endpoint.

## 2. The basepoint covariance lemma

### Lemma 2.1 (a median controls mean by standard deviation)

Let \(T\) be any square-integrable real random variable for which zero is a
median.  If \(m=\mathbb ET\) and
\(\sigma^2=\operatorname {Var}T\), then

\[
                               \boxed{|m|\le\sigma}.    \tag{2.1}
\]

#### Proof

Suppose first that \(m>0\).  Cantelli's one-sided Chebyshev inequality
gives

\[
 \mathbb P(T\le0)
 =\mathbb P(T-m\le-m)
 \le {\sigma^2\over \sigma^2+m^2}.                     \tag{2.2}
\]

The left side is at least one half by the median hypothesis.  Hence
\(m^2\le\sigma^2\).  Apply the same argument to \(-T\) when \(m<0\).
The case \(m=0\) is immediate. \(\square\)

### Theorem 2.2 (basepoint covariance)

Under (1.1)--(1.3),

\[
 \boxed{
 \begin{aligned}
  \mathbb E\,\sigma_q^2N_qN_q^T&\preceq I,\\
  \operatorname {Cov}(B_q)&\preceq I,\\
  \operatorname {Cov}(Y_q)&\preceq4I,\\
  \mathbb E Y_qY_q^T&\preceq4I.
 \end{aligned}}                                             \tag{2.3}
\]

#### Proof

The law of total covariance is

\[
 I=\operatorname {Cov}(X)
  =\mathbb E\operatorname {Cov}(X\mid q)
     +\operatorname {Cov}(\mathbb E[X\mid q])
  =\mathbb E\sigma_q^2N_qN_q^T+\operatorname {Cov}(B_q).
                                                               \tag{2.4}
\]

Both summands are positive semidefinite, proving the first two lines.
Write \(Z_q=m_qN_q\), so that \(Y_q=B_q-Z_q\).  For every vector \(v\),

\[
\begin{aligned}
 \operatorname {Var}\langle v,Y_q\rangle
 &\le2\operatorname {Var}\langle v,B_q\rangle
       +2\operatorname {Var}\langle v,Z_q\rangle\\
 &\le2|v|^2+2\mathbb E[m_q^2\langle v,N_q\rangle^2]\\
 &\le2|v|^2+2\mathbb E[\sigma_q^2\langle v,N_q\rangle^2]
 \le4|v|^2,                                               \tag{2.5}
\end{aligned}
\]

where Lemma 2.1 and (2.4) were used in the last line.  This proves the
third assertion.

There is a slightly sharper uncentered estimate.  Since
\(\mathbb EB_q=\mathbb EX=0\), (2.4) gives
\(\mathbb E B_qB_q^T\preceq I\).  Also, with \(Z_q=m_qN_q\),
Lemma 2.1 and (2.4) give \(\mathbb E Z_qZ_q^T\preceq I\).
The pointwise matrix inequality

\[
 (B_q-Z_q)(B_q-Z_q)^T
 \preceq2B_qB_q^T+2Z_qZ_q^T
\]

therefore proves
\(\mathbb E Y_qY_q^T\preceq4I\).  The covariance bound follows
immediately. \(\square\)

The constant four is all that the proposed projection argument needs.
The important limitation is that (2.3) is an operator bound, not a trace
bound: it permits \(\mathbb E|Y|^2\) to be of order \(n\).

## 3. A fixed quotient mass lies at the witness scale

Use the median version

\[
 D_1(\mu)=\sup_{\operatorname {Lip}(g)\le1}
       \int|g-\operatorname {med}g|\,d\mu.              \tag{3.1}
\]

For log-concave probabilities on Euclidean space, E. Milman's equivalence
theorem applies to locally Lipschitz functions and gives a universal
constant \(C_M\) such that

\[
                         C_P(\mu)\le C_MD_1(\mu)^2.     \tag{3.2}
\]

The theorem applies on the affine support and is invariant under translation
and Euclidean isometry.  In the present full-dimensional isotropic setting
there is no degeneracy to remove.

### Theorem 3.1 (single-scale quotient packet)

Let \(f\) satisfy (1.4), set

\[
                         D=\mathbb E|f|=\mathbb E|T_q|,
                         \qquad K=D_1(\mu),              \tag{3.3}
\]

and suppose \(D\ge\kappa K\), where \(0<\kappa\le1\).  Then the quotient
set

\[
 G=\left\{q:{D\over2\sqrt2}\le\sigma_q
                \le {4C_M\over\kappa^2}D\right\}       \tag{3.4}
\]

satisfies

\[
                              \boxed{\eta(G)\ge
                                  {\kappa^2\over16C_M}}. \tag{3.5}
\]

#### Proof

Lemma 2.1 gives, conditionally on \(q\),

\[
 \mathbb E_q|T_q|
 \le (\mathbb E_qT_q^2)^{1/2}
 = (\sigma_q^2+m_q^2)^{1/2}
 \le\sqrt2\,\sigma_q.                                  \tag{3.6}
\]

Consequently

\[
                         \mathbb E\sigma_q\ge {D\over\sqrt2}
                         \ge {\kappa K\over\sqrt2}.     \tag{3.7}
\]

The weak gradient of a signed-distance function has norm at most one.
It is in \(L^2(\mu)\) because \(\mu\) has a second moment.  Applying
(3.2) and then total variance gives

\[
 \mathbb E\sigma_q^2
 \le\operatorname {Var}(T_q)
 =\operatorname {Var}_\mu f
 \le C_MK^2.                                            \tag{3.8}
\]

Paley--Zygmund applied to the nonnegative variable \(\sigma_q\) yields

\[
 \mathbb P\left\{\sigma_q\ge\frac12\mathbb E\sigma_q\right\}
 \ge { (\mathbb E\sigma_q)^2\over4\mathbb E\sigma_q^2}
 \ge {\kappa^2\over8C_M}.                              \tag{3.9}
\]

On this event \(\sigma_q\ge D/(2\sqrt2)\).  Since
\(K\le D/\kappa\), Markov's inequality gives

\[
 \mathbb P\left\{\sigma_q>{4C_M\over\kappa^2}D\right\}
 \le {C_MK^2\over(4C_MD/\kappa^2)^2}
 \le {\kappa^2\over16C_M}.                             \tag{3.10}
\]

Subtracting (3.10) from (3.9) proves (3.5). \(\square\)

### Lemma 3.2 (two bounded endpoint slabs)

There are universal constants \(c_0,C_0,\beta_0>0\) such that, for every
one-dimensional log-concave probability with median zero and standard
deviation \(\sigma\),

\[
\begin{aligned}
 \mathbb P\{c_0\sigma\le T\le C_0\sigma\}&\ge\beta_0,\\
 \mathbb P\{-C_0\sigma\le T\le-c_0\sigma\}&\ge\beta_0.
                                                               \tag{3.11}
\end{aligned}
\]

#### Proof

Use the standard one-dimensional log-concave density estimate
\(\|q\|_\infty\le C_d/\sigma\).  Take
\(c_0=(8C_d)^{-1}\) and \(C_0=4\).  The mass between zero and
\(c_0\sigma\) is at most \(1/8\).  By Lemma 2.1,
\(\mathbb ET^2\le2\sigma^2\), so
\(\mathbb P(|T|>4\sigma)\le1/8\).  Each side of zero has mass one half;
therefore each interval in (3.11) has mass at least \(1/4\).  One may take
\(\beta_0=1/4\). \(\square\)

Combining Theorem 3.1 and Lemma 3.2 produces fixed positive bulk mass in
two endpoint slabs whose distances from the zero interface lie between
universal multiples of \(D\) (with constants depending only on the fixed
near-maximality parameter \(\kappa\)).

## 4. The orientation tensor gives dispersion, not capture

Assume now exact conditional balance, so
\(\mathbb E_q\operatorname {sgn}(T_q)=0\).  Put

\[
                       h(X)=\operatorname {sgn}(T_q),
                       \qquad a_q=\mathbb E_q|T_q|.     \tag{4.1}
\]

Then

\[
\begin{aligned}
 \mathbb E[hX\otimes N\mid q]
 &=Y_q\otimes N_q\,\mathbb E_qh
     +N_q\otimes N_q\,\mathbb E_q(hT_q)\\
 &=a_qN_q\otimes N_q.                                  \tag{4.2}
\end{aligned}
\]

Thus (0.6) also equals \(\mathbb E[hX\otimes N]\).

### Lemma 4.1 (transport tensor bound)

The matrix \(M\) is positive semidefinite and

\[
                         \operatorname {tr}M=D,
                         \qquad \|M\|_{HS}\le1.        \tag{4.3}
\]

#### Proof

Positivity and the trace identity follow from (0.6).  For every matrix
\(A\), isotropy and \(|h|=|N|=1\) give

\[
\begin{aligned}
 |\langle A,M\rangle_{HS}|
 &=|\mathbb E[hX^TAN]|\\
 &\le\bigl(\mathbb E|A^TX|^2\bigr)^{1/2}
      \bigl(\mathbb E h^2|N|^2\bigr)^{1/2}
 =\|A\|_{HS}.                                           \tag{4.4}
\end{aligned}
\]

Hilbert--Schmidt duality proves the norm bound. \(\square\)

In particular,

\[
 { (\operatorname {tr}M)^2\over\operatorname {tr}(M^2)}
 \ge D^2.                                                \tag{4.5}
\]

The same conclusion is carried by the scale packet.  Lemma 3.2 gives
\(a_q\ge c\sigma_q\) on \(G\), so

\[
 \operatorname {tr}M_G
 :=\operatorname {tr}\int_Ga_qN_qN_q^T\,d\eta(q)
 \ge c_\kappa D,                                        \tag{4.6}
\]

while \(0\preceq M_G\preceq M\) and
\(\|M_G\|_{HS}\le\|M\|_{HS}\le1\).  Hence the normalized good-packet
direction matrix has effective rank at least \(c_\kappa D^2\).

### Proposition 4.2 (the low-dimensional capture inference is false)

There is no function of (4.3), (4.5), and the basepoint bound (2.3) which
produces a rank-\(CD^2\) projection \(P\) satisfying

\[
                       \operatorname {tr}(PM)\ge cD    \tag{4.7}
\]

with universal \(c,C>0\).

#### Proof

For any integer \(N\ge D^2\), the matrix

\[
                              M_N={D\over N}I_N         \tag{4.8}
\]

has trace \(D\) and Hilbert--Schmidt norm \(D/\sqrt N\le1\).  For every
rank-\(r\) projection,

\[
                  {\operatorname {tr}(PM_N)\over
                         \operatorname {tr}M_N}={r\over N}. \tag{4.9}
\]

Taking \(N/(CD^2)\to\infty\) contradicts (4.7).  The basepoint covariance
bound is compatible with this example: a spherical basepoint law of radius
at most \(2\sqrt N\) has covariance at most \(4I\). \(\square\)

A spectral threshold does not repair the direction of the implication.  If
\(H_\tau\) is spanned by the eigenvectors of \(M\) with eigenvalue at least
\(\tau/D\), then

\[
                         \dim H_\tau\le {D^2\over\tau^2}. \tag{4.10}
\]

But (4.8) has \(H_\tau=\{0\}\) whenever \(N>D^2/\tau\), even though its
entire trace is present.  Thus (4.10) bounds the dimension of a chosen
spectral part; it does not prove that the part carries any fixed trace.

There is also a contradiction-framework obstruction to deriving capture
from constant-factor near maximality alone.  Suppose hypothetically that an
isotropic log-concave block \(\nu\) has a centered 1-Lipschitz witness
\(g\) with first-moment scale \(D\), variance comparable with \(D^2\), and
nonlinear regression residual

\[
 \mathbb E\left[g(X)-\langle \mathbb E(Xg),X\rangle\right]^2
 \ge\theta^2D^2.                                        \tag{4.11}
\]

For independent blocks set

\[
 F_m(X_1,\ldots,X_m)=m^{-1/2}\sum_{j=1}^mg(X_j).        \tag{4.12}
\]

Then \(\nu^{\otimes m}\) is isotropic and log-concave, \(F_m\) is
1-Lipschitz, and the central limit theorem gives
\(\mathbb E|F_m|\asymp D\).  Bobkov--Houdré tensorization together with
Milman's equivalence keeps \(D_1(\nu^{\otimes m})\asymp D_1(\nu)\), so a
constant-factor near-maximal block produces a constant-factor near-maximal
tensor witness.

For every fixed \(k\), a joint triangular-array central limit argument for
\((F_m,P(X_1,\ldots,X_m))\), uniform over rank-\(k\) projections after
separating their finitely many high-leverage blocks, gives

\[
 \liminf_{m\to\infty}
 \inf_{\substack{\operatorname {rank}P\le k\\
                  \operatorname {Lip}\phi\le1}}
 \mathbb E|F_m-\phi(PX)|
 \ge\sqrt{2\over\pi}\,\theta D.                        \tag{4.13}
\]

In the Gaussian limit, (4.13) is simply the least-absolute-deviation error
of a conditional Gaussian whose unexplained variance is at least the
residual in (4.11).  Taking \(k=CD^2\) and then \(m\) large shows that
near maximality, isotropy, and tensorization do not imply a rank-\(O(D^2)\)
ridge description.

After translating \(F_m\) by one of its medians, replacing it by the signed
distance to its zero interface increases its pointwise absolute value and
preserves 1-Lipschitzness, so it also preserves the median first-moment
scale.  What this replacement does **not** prove
is cellwise balance or stationarity of the new interface.  Consequently a
valid capture theorem would have to use the exact Euler equation and global
maximality in a way which excludes this tensorized nonlinear residual.  It
cannot follow from the covariance and transport estimates established
above.

## 5. What a genuinely captured packet would imply

Although capture is not proved, it is useful to state exactly what it would
buy.  Fix constants \(\alpha,a,b,\varepsilon>0\).  Suppose \(G_1\subset Q\)
has \(\eta(G_1)\ge\alpha\), and for every \(q\in G_1\)

\[
\begin{gathered}
 \mathbb P_q\{aD\le T_q\le bD\}\ge\alpha,
 \qquad
 \mathbb P_q\{-bD\le T_q\le-aD\}\ge\alpha,            \tag{5.1}\\
                         |P_{H^\perp}N_q|\le\varepsilon,\qquad
                         d:=\dim H\le C_HD^2.            \tag{5.2}
\end{gathered}
\]

Theorem 3.1 and Lemma 3.2 provide (5.1); (5.2) is the unproved capture
input.  Let \(\mu_+\) and \(\mu_-\) be the restrictions of \(\mu\) to the
two slab packets in (5.1), trimmed if necessary to equal mass
\(\rho\ge\alpha^2\).  Write

\[
                         \nu_\pm=(P_H)_\#\mu_\pm.       \tag{5.3}
\]

### 5.1 A precise projected matching dichotomy

For finite measures \(\nu_+,\nu_-\) of equal mass and \(L>0\), define
\(\mathfrak C_L(\nu_+,\nu_-)\) to be the largest mass of a partial
coupling with marginals dominated by \(\nu_+,\nu_-\) and supported on

\[
                         R_L=\{(z,z'):|z-z'|\le L\}.     \tag{5.4}
\]

The weighted bipartite matching theorem, first for finite atomic
approximations and then by weak compactness, gives the following alternative.

* Either \(\mathfrak C_L\ge\rho/2\);
* or there are submeasures \(\bar\nu_\pm\le\nu_\pm\), each of mass at
  least \(\rho/2\), whose essential supports have mutual distance at least
  \(L\).

For completeness, in the finite atomic case maximum matching equals minimum
weighted vertex cover.  If the matching mass is below \(\rho/2\), remove a
cover of total weight below \(\rho/2\).  The remaining mass on each side is
at least \(\rho/2\), and no remaining cross edge belongs to \(R_L\).
Approximate general finite measures on increasing compact sets by finite
partitions, use \(R_{L+\delta}\), pass to a weakly convergent subsequence,
and let \(\delta\downarrow0\).  This proves the stated essential-support
version without an unproved measurable-selection assertion.

### Theorem 5.1 (the separated branch closes dimension-freely)

If the second alternative holds with \(L=cD\), then \(D\) is bounded by a
constant depending only on \(c,\rho,C_H\) and the universal constant in the
known \(\sqrt{\log d}\) KLS theorem.

#### Proof

Let \(A_+,A_-\subset H\) be essential supports of the two residual
submeasures.  Define

\[
             \phi(z)={d(z,A_-)-d(z,A_+)\over2}.          \tag{5.5}
\]

The function is 1-Lipschitz, is at least \(L/2\) on \(A_+\), and at most
\(-L/2\) on \(A_-\).  If \(Z,Z'\) are independent with law
\(\bar\mu=(P_H)_\#\mu\), the two cross events have total probability at
least \(2(\rho/2)^2\).  Therefore

\[
 \mathbb E|\phi(Z)-\mathbb E\phi(Z)|
 \ge\frac12\mathbb E|\phi(Z)-\phi(Z')|
 \ge {\rho^2L\over4}.                                  \tag{5.6}
\]

The marginal \(\bar\mu\) is isotropic and log-concave on \(H\).  The known
dimension-dependent KLS theorem gives

\[
 \psi_{\bar\mu}\ge {c_K\over\sqrt{\log(e+d)}}.
                                                               \tag{5.7}
\]

Cheeger's inequality then gives

\[
 C_P(\bar\mu)\le {4\over c_K^2}\log(e+d),
 \qquad
 \mathbb E|\phi-\mathbb E\phi|
 \le {2\over c_K}\sqrt{\log(e+d)}.                    \tag{5.8}
\]

Combining (5.6), \(L=cD\), and \(d\le C_HD^2\),

\[
 D\le {8\over c_Kc\rho^2}
        \sqrt{\log(e+C_HD^2)}.                          \tag{5.9}
\]

This scalar self-consistency inequality bounds \(D\) universally.  For
example, square it, use
\(e+C_HD^2\le(e+C_H)(1+D^2)\) and
\(\log(1+D^2)\le D\), and solve the resulting quadratic inequality in
\(D\). \(\square\)

This is stronger and cleaner than a Paouris classification: once a
low-dimensional projected witness survives, the already known
\(\sqrt{\log d}\) theorem closes it.

## 6. Quantitative projected collision creates transverse rank

Assume now that the first matching alternative holds.  Lift a partial
coupling of \(\nu_+,\nu_-\) to a coupling \(\pi\) of the original slab
submeasures by using regular conditional laws over \(P_HX\).  Its mass
\(\gamma\) is at least \(\rho/2\), and

\[
                         |P_H(X_+-X_-)|\le L           \tag{6.1}
\]

for \(\pi\)-almost every pair.  Since \(f(X_+)\ge aD\),
\(f(X_-)\le-aD\), and \(f\) is 1-Lipschitz,

\[
 |X_+-X_-|\ge2aD.                                       \tag{6.2}
\]

Choose \(L=aD\).  Equations (6.1)--(6.2) imply

\[
                  |P_{H^\perp}(X_+-X_-)|\ge\sqrt3aD.   \tag{6.3}
\]

Write \(X_+=Y_q+tN_q\) and \(X_-=Y_{q'}+t'N_{q'}\), where
\(|t|,|t'|\le bD\).  Under (5.2),

\[
 |P_{H^\perp}(Y_q-Y_{q'})|
 \ge(\sqrt3a-2b\varepsilon)D.                          \tag{6.4}
\]

Take \(\varepsilon\le(\sqrt3-1)a/(2b)\), and set

\[
 U(q,q')={P_{H^\perp}(Y_q-Y_{q'})
              \over|P_{H^\perp}(Y_q-Y_{q'})|},
 \qquad
 Q_c=\int U\otimes U\,d\pi.                            \tag{6.5}
\]

### Theorem 6.1 (collision-to-rank, not collision-to-saving)

The collision matrix satisfies

\[
 \boxed{
     \operatorname {tr}Q_c=\gamma,
     \qquad
     \|Q_c\|_{op}\le {20\over a^2D^2},
     \qquad
     {\operatorname {tr}Q_c\over\|Q_c\|_{op}}
       \ge {\gamma a^2D^2\over20}.}                    \tag{6.6}
\]

#### Proof

The quotient marginals of the lifted coupling are dominated by \(\eta\):
the slab-selection probabilities are at most one.  Theorem 2.2 therefore
gives, for every unit vector \(v\),

\[
\begin{aligned}
 a^2D^2\langle v,Q_cv\rangle
 &\le\int\langle v,Y_q-Y_{q'}\rangle^2\,d\pi\\
 &\le2\int\langle v,Y_q\rangle^2\,d\pi_1
       +2\int\langle v,Y_{q'}\rangle^2\,d\pi_2\\
 &\le4\,\mathbb E_\eta\langle v,Y_q\rangle^2
 \le20.                                                  \tag{6.7}
\end{aligned}
\]

The first equality in (6.6) follows from \(|U|=1\); the other two follow
from (6.7). \(\square\)

The conclusion is compatible with every covariance constraint.  It says
that collision in the projection is paid for by a new collection of
physical displacements spread through at least \(cD^2\) transverse
directions.  It does **not** say that two interface patches are physically
close.  In fact (6.3) says the reverse.

Any local chord replacement, bevel, or tube-splicing argument needs close
points in the original Euclidean metric.  A projected matching supplies
only close images.  To turn the latter into a competitor, one would have to
solve an isoperimetric gluing problem in the \(H^\perp\)-fibers and control
the transverse graph perimeter and volume repair.  A uniform solution of
that problem is a marginal-to-whole KLS transfer; it is not a consequence
of (0.2).

One can append the span of the new displacement packet to \(H\), but
Theorem 6.1 is again a lower-effective-rank statement.  In an ambient space
of dimension \(n\), the covariance trace budget permits order \(n/D^2\)
orthogonal packets.  Thus this iteration has dimension-dependent depth and
cannot be used in a dimension-free proof.

## 7. Why Paouris and thin shell do not classify the collision branch

Suppose optimistically that a surviving marginal has dimension
\(d\asymp D^2\).  The separated tail packets in Theorem 5.1 are then at
distance comparable with \(\sqrt d\).  Thin shell places a fixed fraction
of an isotropic log-concave marginal in an annulus of radius \(\sqrt d\)
and universal width.  After discarding a small fixed amount of each tail,
a \(c\sqrt d\) separation is therefore predominantly angular, not radial.

This observation excludes the two naive classifications rather than proving
one of them.

* A fixed-mass radial quantile gap is \(O(1)\) by the thin-shell variance
  bound for \(|Z|\).
* A fixed-mass affine quantile gap is \(O(1)\), because every linear
  functional of an isotropic log-concave vector is a one-dimensional
  log-concave variable of variance one.

Thus a hypothetical \(c\sqrt d\) witness must be nonlinear and angular.
Paouris controls radial tails and thin shell controls radial fluctuation;
neither controls angular checkerboards.  The following elementary
barycenter estimate makes the same point.  If \(A\) has probability
\(\alpha\) under an isotropic law and \(b_A=\mathbb E[X\mid A]\), then

\[
                              |b_A|\le\sqrt{1-\alpha\over\alpha}. \tag{7.1}
\]

Indeed the complement barycenter is
\(-\alpha b_A/(1-\alpha)\), and the variance in direction \(b_A\) is at
least \(\alpha|b_A|^2/(1-\alpha)\).  Far constant-mass packets can still
obey (7.1) by interlacing angularly.  The atomic simplex in Section 8.3 is
an exact shell coloring of this kind.

The selected basepoint law is not log-concave, so Paouris cannot even be
applied to \(Y_q\) or to the collision direction measure in (6.5).  It can
be applied only to the whole projected marginal, which does not retain the
cell labels needed to distinguish the two phases.  A theorem converting
whole-marginal radial control into angular regularity of a maximizing
signed-distance partition would be a new isoperimetric inverse of
conjecture strength.

## 8. Exact model audit

### 8.1 Isotropic cube: the affine survivor

Let

\[
                   \mu=\operatorname {Unif}[-\sqrt3,\sqrt3]^n,
                   \qquad f(x)=x_1.                     \tag{8.1}
\]

The normal cells are

\[
 Y=(0,X_2,\ldots,X_n),\qquad N=e_1,qquad
 T=X_1.                                                  \tag{8.2}
\]

They are balanced, \(m=0\), \(\sigma=1\), and

\[
 D=\mathbb E|T|={\sqrt3\over2},qquad
 \operatorname {Cov}(Y)=I-e_1e_1^T,qquad
 M={\sqrt3\over2}e_1e_1^T.                             \tag{8.3}
\]

The one-dimensional subspace \(H=\operatorname {span}(e_1)\) captures all
directions, the projected phases remain separated, and the projected
witness is exactly \(z\mapsto z\).  This is the first branch of Theorem
5.1 with no loss.

### 8.2 Isotropic radial exponential: diffuse projected collision

Let

\[
 d\mu(x)=c_{n,\lambda}e^{-\lambda|x|}\,dx,
 \qquad \lambda=\sqrt{n+1}.                             \tag{8.4}
\]

Then \(R=|X|\) has the Gamma\((n,\lambda)\) law and

\[
 \mathbb ER={n\over\lambda},\qquad
 \operatorname {Var}R={n\over n+1},                    \tag{8.5}
\]

which makes \(\mu\) isotropic.  Let \(r_0\) be the median of \(R\) and
take \(f(x)=|x|-r_0\).  With \(U=X/|X|\),

\[
                         Y=r_0U,qquad N=U,qquad T=R-r_0. \tag{8.6}
\]

The direction \(U\) is uniform on the sphere and independent of \(T\).
The cells are exactly balanced.  Moreover

\[
 \sigma^2={n\over n+1},qquad |m|\le\sigma,qquad
 c\le D=\mathbb E|R-r_0|\le\sqrt2,                     \tag{8.7}
\]

where the lower bound follows from the universal upper density bound for a
variance-one log-concave law.  The matrices are

\[
 \operatorname {Cov}(Y)={r_0^2\over n}I,qquad
                         M={D\over n}I.                 \tag{8.8}
\]

Thus a rank-\(d\) subspace captures exactly \(d/n\) of the orientation
trace.  For \(d\ll n\), the projections of inner and outer radial slabs
have large overlap: the omitted radius can be varied while the first
\(d\) coordinates are fixed.  The matched full points remain separated in
\(H^\perp\), exactly as in (6.3).

The median sphere is not being asserted to be the global T3 maximizer, but
its signed-distance second variation is semistable; the degree-one
spherical harmonics are exact equality modes.  Hence projected overlap by
itself cannot imply a local negative second variation or an automatic
first-moment saving.  This model realizes the diffuse matrix obstruction
inside an isotropic log-concave law.

### 8.3 Atomic regular simplex: sharp maximal-scale shell coloring

Let \(N=k+1\) be even, and let \(v_1,\ldots,v_N\) be the vertices of a
centered isotropic regular simplex:

\[
 |v_i|^2=k,qquad \langle v_i,v_j\rangle=-1\ (i\ne j),
 \qquad {1\over N}\sum_iv_iv_i^T=I.                    \tag{8.9}
\]

Partition the vertices into two equal classes and match them in pairs
\((p_j,q_j)\).  Every cross distance is
\(L=\sqrt{2N}\).  The function taking values \(L/2\) on the positive
class and \(-L/2\) on the negative class is 1-Lipschitz and globally
maximizes median absolute deviation.  Its signed-distance value is

\[
                              D={L\over2}=\sqrt{N\over2}. \tag{8.10}
\]

Disintegrate each matched pair as

\[
 Y_j={p_j+q_j\over2},qquad
 N_j={p_j-q_j\over L},qquad
 T_j\in\{-D,D\}\text{ equiprobably}.                   \tag{8.11}
\]

The cells are exactly balanced and \(\sigma_j=D\).  Difference vectors
from disjoint matched pairs are orthogonal, because

\[
 \langle p_i-q_i,p_j-q_j\rangle=(-1)-(-1)-(-1)+(-1)=0. \tag{8.12}
\]

There are \(N/2=D^2\) such directions, and

\[
 M={2D\over N}\sum_{j=1}^{N/2}N_jN_j^T,qquad
 \operatorname {tr}M=D,qquad \|M\|_{HS}=1.            \tag{8.13}
\]

This exactly saturates the \(D^2\) rank scale.  Taking \(H\) to be their
span gives a full surviving nonlinear witness in dimension \(D^2\).
All points lie on one thin shell.  Paouris-type radial information and
thin-shell information therefore cannot rule out the coloring.

The failed hypothesis is log-concavity.  This example shows why the last
step cannot be a metric, covariance, or radial classification theorem.

### 8.4 Uniform simplex: convex filling creates contact

For a regular simplex written as

\[
 K=\{(y,z):0\le z\le H,\ y\in(1-z/H)B\},              \tag{8.14}
\]

the half-volume cap \(E=\{z>c\}\) has

\[
                         c=H(1-2^{-1/k}).               \tag{8.15}
\]

On a vertical normal line based at \(y\), the support ends at
\(z_+(y)=H(1-\|y\|_B)\).  The two cell masses, per unit base area, are

\[
                         a(y)=\rho(z_+(y)-c),qquad
                         b(y)=\rho c.                  \tag{8.16}
\]

They are unequal except on one level set.  For basepoints outside the cap
cross-section, the cell has only lower-side mass.  Thus the obvious
uniform-simplex analogue of the atomic coloring fails the Euler balance
equation on a fixed set of normal/contact cells.  Convex filling is not a
small perturbation of the atomic shell model; it supplies exactly the
global log-concave compatibility absent in Section 8.3.

### 8.5 Product exponentials and the maximum interface

Let \(Y_i\) be independent rate-one exponentials and \(X_i=Y_i-1\).  This
law is isotropic and log-concave.  Choose \(z_n\) by

\[
                         (1-e^{-z_n})^n={1\over2},
 \qquad E=\{\max_iY_i\le z_n\}.                         \tag{8.17}
\]

The signed distance from the boundary equals
\(z_n-\max_iY_i\) inside \(E\), while outside it is
\((\sum_i(Y_i-z_n)_+^2)^{1/2}\).  Its expectation is between two universal
positive constants.  The upper bound follows from

\[
 \mathbb E\sum_i(Y_i-z_n)_+^2=2ne^{-z_n}=O(1)           \tag{8.18}
\]

outside, and from
\((1-e^{-s})^n\le e^{-ne^{-s}}\) in the layer-cake integral inside.
A fixed interval next to the median level gives the lower bound.

Permutation symmetry makes the regular-facet orientation matrix
proportional to \(I_n\), so any fixed-dimensional coordinate projection
captures a vanishing fraction.  The omitted coordinates determine the
maximum with fixed probability, and hence the projected phases have fixed
overlap.

This is not a stationary counterexample.  On a regular facet the exact
normal-cell density is

\[
                  w(t)=Ae^{-t}\mathbf1_{[-d,\infty)}(t),
                  \qquad d=z_n-\max_{j\ne i}Y_j,         \tag{8.19}
\]

so the negative and positive masses have ratio
\((e^d-1):1\) and are equal only when \(d=\log2\).  Moreover the mass of
points assigned to ridge normal cones tends to

\[
                              {1-\log2\over2}>0.         \tag{8.20}
\]

A simultaneous bevel of pairwise ridges gives a fixed perimeter saving.
Thus the model shows both failure modes of a projection-only argument:
diffuse projected collision is real, while the actual saving is carried by
full-space ridge geometry which the projection has discarded.

## 9. Exact remaining statement

The proved implication is

\[
\begin{gathered}
 \text{near-maximal balanced signed distance}
 \Longrightarrow
 \operatorname {Cov}(Y)\preceq4I
 \text{ and fixed mass of }\sigma\asymp D,\\
 \Longrightarrow
 \text{orientation effective rank }\Omega(D^2).        \tag{9.1}
\end{gathered}
\]

What is **not** proved is

\[
 \text{effective rank }\Omega(D^2)
 \Longrightarrow
 \text{fixed direction mass in one rank-}O(D^2)\text{ subspace}. \tag{9.2}
\]

Equation (9.2) is false for the available matrix data.  Conditional on a
capturing subspace, the separated projected branch closes by Theorem 5.1.
The collision branch gives only Theorem 6.1.  Therefore the genuinely new
lemma required by this route is the following.

> **Extremal transverse-contact lemma.**  For a globally maximizing
> half-mass signed-distance interface of an isotropic log-concave law, a
> fixed-mass coupling of opposite \(D\)-scale normal slabs whose projection
> collides, and whose basepoint displacement law has effective rank
> \(\Omega(D^2)\), forces a full-space competitor with a universal relative
> first-moment or perimeter saving, unless \(D=O(1)\).

The radial exponential shows that collision alone is insufficient; the
atomic simplex shows that global maximality, isotropy, balance, and thin
shell are insufficient without log-concavity; the product maximum shows
that the saving, when present, is a simultaneous ridge/contact operation.
Proving the displayed lemma would supply exactly the missing global
log-concave compatibility.  It is not a consequence of Milman's theorem,
Paouris, thin shell, the law of total covariance, or the transport tensor,
and assuming it would simply rename the unresolved high-rank branch.
