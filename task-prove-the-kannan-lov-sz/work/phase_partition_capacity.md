# Phase partitions, tilt-law capacity, and the additive-scale obstruction

## Executive verdict

The natural tilt law does give a noncircular, dimension-free capacity
estimate for **separated phase packets**.  At precision

\[
                         t={\tau\over K},
\]

its Poincare constant is at most \(t(1+\tau)\), and hence its Cheeger
constant is at least \(c/\sqrt{t(1+\tau)}\).  Meanwhile the posterior
centroid field is globally \(t^{-1}\)-Lipschitz.  Therefore two central,
near-saturated phase packets of masses at least \(a\), whose unoriented
directions are separated by \(\gamma\), force an intervening tilt-law mass

\[
                         \gtrsim_\delta
 a\min\left\{1,{\gamma\over\sqrt{1+\tau}}\right\}.       \tag{0.1}
\]

If every intervening state is central and has centroid defect at least
\(e\), (0.1) is charged to the localization-profile derivative.
This proves a rigorous phase-cell capacity theorem and works well for
product/facet phase cells.

It does **not** turn high Fisher effective rank into a KLS closure.  There
are three independent obstructions.

1. Effective rank alone need not give two fixed-mass, angularly separated
   packets.  The uniform projective-spherical phase law has effective rank
   \(m\), but spherical isoperimetry puts every two sets of fixed positive
   mass within angular distance \(O(m^{-1/2})\).  The radial Gaussian test
   realizes this behavior exactly.
2. Conditioning the tilt law on central, low-defect states destroys
   log-concavity and need not preserve any Poincare or Cheeger bound.
   Applying expansion to that weighted phase law would insert an unproved
   KLS statement.
3. Even under the optimistic assumption that every diffuse phase family
   is converted into a Gaussian factor, the audited angular modulus requires

\[
 e_\tau\lesssim_\delta
 {(1+\tau)^{-6}\over\sqrt{\log(e+\tau)}}.               \tag{0.2}
\]

The exact profile budget is logarithmic:

\[
 \int q(\tau){d\tau\over\tau}=O(1),\qquad
 q={\mathbb E[I(g)(1-\eta^2)]\over\mathbb E I(g)}.       \tag{0.3}
\]

Every polynomial threshold in (0.2), including the extra
\((1+\tau)^{-1/2}\) capacity charge, is integrable against
\(d\tau/\tau\).  Thus neither dyadic nor additive \(\tau\)-pigeonholing
forces a scale meeting it.  An explicit admissible scalar profile can stay
above the required threshold at every scale while consuming only finite
total profile growth.

Gaussian parity, radial cuts, and product phase cells realize the three
branches sharply.  Parity shows that two phases may connect through a
noncentral tail corridor whose \(I(g)\)-charge is exponentially small; the
posterior is a full Gaussian factor.  Radial phases saturate the Poincare
and angular threshold without any separator; again they are the Gaussian
factor branch.  Product cells have a large tie separator and are correctly
detected by (0.1), but the scale-selection mismatch remains.

Accordingly, the proposed phase-partition closure is **refuted as a
consequence of high effective rank, the pointwise angular estimate, and the
known profile budget alone**.  A complete dichotomy would additionally need
a new quantitative inverse theorem saying that every diffuse low-defect
phase law comes from one common Gaussian factor, plus a profile mechanism
stronger than logarithmic scale pigeonholing.

## 1. The natural tilt law and its exact scale

Let \(\mu\) have Poincare constant \(K<\infty\), let

\[
 C_t=tX+\sqrt t\,G,\qquad X\sim\mu,\quad G\sim N(0,I),
                                                               \tag{1.1}
\]

and denote the law of \(C_t\) by \(\nu_t\).  For every locally Lipschitz
\(\phi\), product Poincare gives

\[
 \begin{aligned}
 \operatorname {Var}\phi(C_t)
 &\le Kt^2\mathbb E|\nabla\phi(C_t)|^2
      +t\mathbb E|\nabla\phi(C_t)|^2\\
 &=t(1+tK)\mathbb E|\nabla\phi(C_t)|^2.                \tag{1.2}
 \end{aligned}
\]

At \(t=\tau/K\),

\[
 \boxed{\quad C_P(\nu_t)\le t(1+\tau).\quad}            \tag{1.3}
\]

The law \(\nu_t\) is log-concave, because it is the affine image of the
convolution of two log-concave laws.  The Buser--Ledoux implication for a
log-concave probability therefore yields

\[
 \boxed{\quad
 h_{\nu_t}\ge {c_0\over\sqrt{t(1+\tau)}}.
 \quad}                                                  \tag{1.4}
\]

No expansion statement for a conditioned or Fisher-weighted tilt law is
being used.

For a fixed set \(S\), write under the posterior \(\pi_c\)

\[
 g(c)=\pi_c(S),\qquad
 v(c)=\operatorname {Cov}_{\pi_c}({\bf1}_S,X),\qquad
 D(c)=\mathbb E_{\pi_c}
 [({\bf1}_S-g)(X-m)(X-m)^T].                            \tag{1.5}
\]

Differentiation of the exponential tilt gives the exact identities

\[
                         \nabla g=v,\qquad\nabla v=D.   \tag{1.6}
\]

The following global bound is useful because it does not require centroid
near-equality.

**Lemma 1.1 (global centroid Lipschitz bound).**  For every posterior state,

\[
                         \|D(c)\|_{\rm op}\le {1\over t}.
                                                               \tag{1.7}
\]

**Proof.**  For a unit vector \(\theta\), put
\(Y=\langle X-m,\theta\rangle\).  The posterior is
\(t\)-strongly log-concave, so its Poincare constant is at most \(1/t\)
and \(\mathbb E Y^2\le1/t\).  Hence

\[
 \operatorname {Var}(Y^2)
 \le {1\over t}\mathbb E|2Y\theta|^2
 \le {4\over t^2}.
\]

Covariance Cauchy--Schwarz and
\(\sqrt{g(1-g)}\le1/2\) give

\[
 |\theta^TD\theta|
 \le\sqrt{g(1-g)}\sqrt{\operatorname {Var}(Y^2)}
 \le {1\over t}.
\]

Taking the supremum of the absolute quadratic form proves (1.7).
\(\square\)

Thus

\[
                         |v(c)-v(c')|\le {|c-c'|\over t}. \tag{1.8}
\]

## 2. Angular calculus on good posterior states

Fix \(\delta\in(0,1/2)\) and \(e\in[0,1/2]\).  On states with
\(v\ne0\), put

\[
 u={v\over|v|},\qquad P_u=uu^T,\qquad
 \eta={\sqrt t|v|\over I(g)},\qquad
 \varepsilon=1-\eta.                                  \tag{2.1}
\]

Define

\[
 G_{\delta,e}
 =\{c:\delta\le g(c)\le1-\delta,\ \varepsilon(c)\le e\}.
                                                               \tag{2.2}
\]

Let

\[
                         i_\delta=
 \min_{r\in[\delta,1-\delta]}I(r)>0.
\]

On \(G_{\delta,e}\),

\[
                         |v|\ge{(1-e)i_\delta\over\sqrt t}.
                                                               \tag{2.3}
\]

Differentiating the normalized vector gives

\[
                         \nabla u={P_uD\over|v|}.       \tag{2.4}
\]

The general angular-stability theorem therefore implies

\[
 \boxed{\quad
 \|\nabla u\|_{HS}^2
 \le {C_\delta\over t}\Omega_\delta(e),\qquad
 \|\nabla P_u\|_{HS}^2
 =2\|\nabla u\|_{HS}^2.
 \quad}                                                  \tag{2.5}
\]

Here

\[
 \Omega_\delta(r)
 \le C_\delta r^{1/6}\{\log(e/r)\}^{1/12}              \tag{2.6}
\]

for all sufficiently small \(r\).

There is also a global separation statement which uses (1.8), rather than
integrating (2.5) along a path that might leave the good set.

**Lemma 2.1 (angular separation forces tilt separation).**  If
\(c,c'\in G_{\delta,e}\) and

\[
                         \|P_{u(c)}-P_{u(c')}\|_{HS}
 \ge\gamma,
                                                               \tag{2.7}
\]

then

\[
 \boxed{\quad
 |c-c'|\ge{(1-e)i_\delta\gamma\over\sqrt2}\sqrt t.
 \quad}                                                  \tag{2.8}
\]

**Proof.**  The component of \(v(c')\) perpendicular to \(u(c)\) has
length

\[
 |v(c')|\sqrt{1-\langle u(c),u(c')\rangle^2}
 ={ |v(c')|\over\sqrt2}
   \|P_{u(c)}-P_{u(c')}\|_{HS}.
\]

Since \(v(c)\) has no such component, (2.3) gives

\[
 |v(c)-v(c')|
 \ge{(1-e)i_\delta\gamma\over\sqrt{2t}}.
\]

Combine this with (1.8). \(\square\)

## 3. A rigorous phase-packet capacity theorem

The next lemma is the positive result of the phase-partition route.

**Theorem 3.1 (capacity of two separated phase packets).**  Let
\(E,F\subset G_{\delta,e}\) be Borel sets satisfying

\[
 \nu_t(E)\ge a,\qquad \nu_t(F)\ge a,\qquad
 \inf_{c\in E,c'\in F}
 \|P_{u(c)}-P_{u(c')}\|_{HS}\ge\gamma.                 \tag{3.1}
\]

Put \(d=\operatorname {dist}(E,F)\),

\[
 r=\min\{d/2,(2h_{\nu_t})^{-1}\},\qquad
 \Sigma=E_r\setminus E.                                \tag{3.2a}
\]

Then \(\Sigma\subset\mathbb R^n\setminus(E\cup F)\) and

\[
 \boxed{\quad
 \nu_t(\Sigma)\ge c_\delta a
       \min\left\{1,{\gamma\over\sqrt{1+\tau}}\right\}.
 \quad}                                                  \tag{3.2}
\]

If, in addition,

\[
 \Sigma\subset\{c:\delta\le g(c)\le1-\delta,\
                    \varepsilon(c)\ge e\},             \tag{3.3}
\]

then the profile dissipation

\[
                         Q(t)=
 \mathbb E_{\nu_t}[I(g)(1-\eta^2)]
                                                               \tag{3.4}
\]

satisfies

\[
 \boxed{\quad
 Q(t)\ge c_\delta\,a e
       \min\left\{1,{\gamma\over\sqrt{1+\tau}}\right\}.
 \quad}                                                  \tag{3.5}
\]

**Proof.**  By Lemma 2.1,

\[
 d:=\operatorname {dist}(E,F)
 \ge c_\delta\gamma\sqrt t.                            \tag{3.6}
\]

For
\(0<s<r\), the \(s\)-neighborhood \(E_s\) is disjoint from \(F\), so

\[
                         a\le\nu_t(E_s)\le1-a.
\]

The defining Cheeger differential inequality gives

\[
 {d\over ds}\nu_t(E_s)\ge h_{\nu_t}a
\]

for almost every \(s<r\).  Therefore

\[
 \nu_t(\Sigma)=\nu_t(E_r)-\nu_t(E)
 \ge h_{\nu_t}ar
 \ge {a\over2}\min\{h_{\nu_t}d,1\}.                   \tag{3.7}
\]

Equations (1.4) and (3.6) prove (3.2).  Under (3.3),
\(I(g)\ge i_\delta\) and
\(1-\eta^2\ge1-\eta=\varepsilon\ge e\) on \(H\).
Thus \(Q(t)\ge i_\delta e\nu_t(\Sigma)\), proving (3.5).
\(\square\)

The assumption (3.3) is substantive.  The complement may instead contain
low-defect intermediate directions, or it may leave the central range
where \(I(g)\) becomes very small.  Neither possibility is charged by
(3.5).

## 4. What global Poincare proves when there are no holes

There is a clean coherence theorem if the phase map satisfies the angular
estimate on the whole tilt space.  It is included to identify the exact
power of the required defect.

Let

\[
                         M=\mathbb E_{\nu_t}P_u.
\]

Since \(\operatorname {tr}P_u=1\),

\[
 \operatorname {Var}_{\nu_t}(P_u)
 :=\mathbb E\|P_u-M\|_{HS}^2
 =1-\operatorname {tr}(M^2).                          \tag{4.1}
\]

**Proposition 4.1 (hole-free phase coherence).**  Suppose
\(\nu_t(G_{\delta,e})=1\).  Then

\[
 \boxed{\quad
 1-\operatorname {tr}(M^2)
 \le C_\delta(1+\tau)\Omega_\delta(e).
 \quad}                                                  \tag{4.2}
\]

Consequently, if the right side is at most \(\alpha<1\), there is a unit
vector \(\theta\) such that

\[
                         \mathbb E\langle u,\theta\rangle^2
 =\lambda_{\max}(M)\ge1-\alpha.                        \tag{4.3}
\]

**Proof.**  Apply the vector-valued Poincare inequality (1.3), component by
component, to \(P_u\), and use (2.5):

\[
 1-\operatorname {tr}(M^2)
 \le t(1+\tau)\mathbb E\|\nabla P_u\|_{HS}^2
 \le C_\delta(1+\tau)\Omega_\delta(e).
\]

Also
\(\lambda_{\max}(M)\ge\operatorname {tr}(M^2)\), because the eigenvalues
are nonnegative and sum to one. \(\square\)

For isotropic \(\mu\), the Fisher phase matrix is

\[
 R_t=\mathbb E_{\nu_t}
 \left[{\eta^2I(g)^2\over g(1-g)}P_u\right]\preceq tI.
                                                               \tag{4.4}
\]

On \(G_{\delta,e}\), its scalar weight is bounded above and below by
positive constants depending only on \((\delta,e)\).  Hence, if a fixed
amount of \(\operatorname {tr}R_t\) lies in \(G_{\delta,e}\), the
unweighted phase law there has effective rank \(\gtrsim1/t\).  If the good
set were the whole space, (4.2) would contradict that rank whenever

\[
                         \Omega_\delta(e)\ll{1\over1+\tau}. \tag{4.5}
\]

For the audited modulus (2.6), a sufficient numerical choice is

\[
 \boxed{\quad
 e=e_\tau:=
 {c_\delta\over(1+\tau)^6\sqrt{\log(e+\tau)}}.
 \quad}                                                  \tag{4.6}
\]

Indeed, after decreasing \(c_\delta\),

\[
 e_\tau^{1/6}
 \{\log(e/e_\tau)\}^{1/12}
 \le {c'_\delta\over1+\tau}.                           \tag{4.7}
\]

The hole-free assumption cannot be replaced by conditioning on
\(G_{\delta,e}\).  That set need not be convex or connected, and the
conditioned law need not be log-concave.  No bound for its Poincare
constant follows from (1.3).  Assuming one would be precisely the forbidden
weighted-phase expansion step.

## 5. Effective rank does not create fixed-mass separated packets

Let \(\sigma_{m-1}\) be uniform measure on the unit sphere and identify
antipodes, so the phase is \(P_u=uu^T\).  Then

\[
                         \mathbb E_{\sigma_{m-1}}P_u={1\over m}I_m,
                                                               \tag{5.1}
\]

which has effective rank \(m\).  Nevertheless, for every fixed
\(a\in(0,1/2)\), spherical isoperimetry gives a constant \(C_a\) such that
any two Borel sets \(A,B\) of projective-spherical measure at least \(a\)
satisfy

\[
                         \operatorname {dist}_{\rm angle}(A,B)
 \le {C_a\over\sqrt m}.                                \tag{5.2}
\]

To see this, the \(r\)-neighborhood of a spherical cap minimizes
neighborhood measure.  Gaussian concentration on the sphere says that an
\(r=C_a/\sqrt m\) neighborhood of any set of mass \(a\) has complement of
mass less than \(a\).  It must therefore meet \(B\).

Thus high rank is compatible with a diffuse continuum of phases and gives
no fixed \((a,\gamma)\) to which Theorem 3.1 applies.  For the KLS scaling
\(m\asymp K/\tau=1/t\), the largest forced angular separation between two
fixed-mass phase sets is only \(O(\sqrt t)\).  Lemma 2.1 then supplies a
tilt separation of order \(t\), while the tilt-law Poincare length is
\(\sqrt{t(1+\tau)}\).  The resulting capacity ratio is only

\[
                         O\left(\sqrt{t\over1+\tau}\right),
                                                               \tag{5.3}
\]

which retains forbidden \(K\)-dependence.

Calling every such diffuse phase law a Gaussian factor is an additional
inverse theorem, not a consequence of effective rank or (2.5).  The radial
test below shows why a Gaussian-factor alternative is necessary, but it
does not prove that it is the only diffuse realization.

Quantitatively, the needed conclusion would have to produce a fixed
subspace \(U\), of dimension at least two (and ultimately comparable to the
phase effective rank).  If \(G\sim N(0,I)\) and \(T_c\) denotes the
Brenier map from \(G\) to the physical posterior, it would have to obey on
a positive set of tilts

\[
 \mathbb E\left|
 T_c(G)-\{m_c+t^{-1/2}P_UG+R_c(P_{U^\perp}G)\}\right|^2
 \le o(t^{-1}).                                         \tag{5.4}
\]

The pointwise splitting theorem supplies (5.4) only with the
one-dimensional space \(U=\operatorname {span}\{u(c)\}\) depending on
\(c\).  Synchronizing those spaces is exactly the missing inverse step.

## 6. Additive \(\tau\)-scales do not repair the profile mismatch

Put

\[
 J(t)=\mathbb E_{\nu_t}I(g),\qquad
 F(t)=\sqrt t\,J(t),\qquad
 q(t)={Q(t)\over J(t)}.                                \tag{6.1}
\]

The exact localization-profile identity is

\[
 {d\log F(t)\over d\log t}={q(t)\over2}.               \tag{6.2}
\]

With \(t=\tau/K\), this becomes

\[
 \boxed{\quad
 \int_{\tau_0}^{\tau_1}q(\tau){d\tau\over\tau}
 =2\log{F(\tau_1/K)\over F(\tau_0/K)}.
 \quad}                                                  \tag{6.3}
\]

On the hypothetical bad-profile plateau, the right side is at most a
numerical constant.  Equation (6.3), not Lebesgue measure \(d\tau\), is the
entire scale budget.

Even the optimistic packet theorem would first choose (4.6).  The
separator charge (3.5) then asks for a relative profile defect smaller than

\[
 q_*(\tau)\asymp_\delta
 {1\over(1+\tau)^{13/2}\sqrt{\log(e+\tau)}}             \tag{6.4}
\]

for fixed packet mass and fixed angular gap.  The exponent \(13/2\) is the
six powers from the angular modulus plus the capacity factor
\((1+\tau)^{-1/2}\).

But

\[
                         \int_1^\infty q_*(\tau)
                         {d\tau\over\tau}<\infty.       \tag{6.5}
\]

Therefore a bounded profile budget cannot force
\(q(\tau)<c q_*(\tau)\) at any scale.  This is not merely a weakness of a
particular pigeonhole partition.  For any small fixed \(b>0\), define the
formal profile

\[
 q_{\rm test}(\tau)=bq_*(\tau),\qquad
 F_{\rm test}(\tau)=F_{\rm test}(1)
 \exp\left\{{1\over2}\int_1^\tau
 q_{\rm test}(s){ds\over s}\right\}.                  \tag{6.6}
\]

Then \(F_{\rm test}\) has a bounded total multiplicative increase, obeys
the exact differential identity (6.2), and remains a fixed multiple of the
required threshold at every \(\tau\).  Replacing dyadic intervals by
additive intervals \([j,j+1]\) cannot distinguish (6.6), because summing
their exact charges reproduces \(d\tau/\tau\).

Even grant the stronger persistence statement that a nonfactor phase
partition survives throughout every unit interval \([j,j+1]\).  The
resulting lower bounds would contribute only

\[
 \sum_{j\ge1}{1\over j}\,j^{-13/2}
 =\sum_{j\ge1}j^{-15/2}<\infty                          \tag{6.7}
\]

to (6.3).  Thus additive windows do not accumulate a divergent charge.
Their apparent linear count is exactly canceled by the \(1/\tau\) weight
in the profile identity.

The same conclusion holds for every positive power
\((1+\tau)^{-p}\).  Thus even a substantially sharper angular modulus with
any finite polynomial inverse would not, by itself, be selected by the
known logarithmic profile budget.  A closing multiscale mechanism must
produce a nonintegrable threshold or an additional additive conserved
quantity.

## 7. Mandatory model tests

### 7.1 Gaussian parity: the noncentral tail corridor

Let the posterior be \(N(m,t^{-1}I_2)\), write its standardized mean as
\((a,b)=\sqrt t\,m\), and take

\[
                         S=\{x_1x_2\ge0\}.
\]

Along the arm \(b=0\),

\[
 g={1\over2},\qquad
 \sqrt t\,v=\varphi(0)(2\Phi(a)-1)e_2,
\]

so

\[
                         u=e_2,\qquad
 \varepsilon=2\Phi(-|a|).                             \tag{7.1}
\]

The orthogonal arm \((0,a)\) has phase \(e_1\) with the same defect.  Their
tilt parameters are separated by

\[
                         \sqrt2\,|a|\sqrt t.           \tag{7.2}
\]

The straight corridor between them passes near standardized mean
\((a/2,a/2)\), where

\[
 g=\Phi(a/2)^2+\Phi(-a/2)^2
 =1-2\Phi(-a/2)\Phi(a/2).                              \tag{7.3}
\]

Hence \(I(g)\le C(1+a)e^{-a^2/8}\).  The two central, almost saturated
phases can therefore communicate through a noncentral corridor carrying
exponentially small profile weight.  Assumption (3.3) fails exactly.

If \(a\asymp\sqrt{1+\tau}\), (7.2) is a fixed fraction of the tilt-law
Poincare length while (7.1) is exponentially smaller than any inverse
power of \(\tau\).  The posterior here is a full Gaussian product, so a
correct structural theorem must put the example in its Gaussian-factor
branch.  Capacity alone cannot charge the corridor.

### 7.2 Radial phases: sharp diffuse saturation

Let

\[
                         \nu=N(0,P I_m),\qquad
 P=t(1+\tau),\qquad u(c)={c\over|c|}.
\]

Then

\[
 \mathbb E uu^T={1\over m}I_m,\qquad
 \|\nabla u(c)\|_{HS}^2={m-1\over|c|^2}.              \tag{7.4}
\]

Since \(|c|^2\asymp mP\) on the Gaussian annulus,

\[
 \mathbb E\|\nabla u\|_{HS}^2\asymp{1\over P}
 ={1\over t(1+\tau)}.                                 \tag{7.5}
\]

The vector Poincare inequality is saturated up to constants:
the phase variance is order one and
\(P\mathbb E\|\nabla u\|^2\asymp1\).  There is no bad separator and no
coherent line.  Comparison with (2.5) shows that the threshold

\[
                         \Omega_\delta(e)\asymp{1\over1+\tau}
                                                               \tag{7.6}
\]

in (4.5) is sharp.  Radial cuts of a Gaussian posterior realize this phase
map; they belong to the full-Gaussian-factor branch and are controlled by
translated thin shell in the original geometric program.

### 7.3 Product phase cells: capacity succeeds but selection fails

Let \(\nu=N(0,P I_m)\), and partition tilt space into the cells

\[
                         \mathcal C_i=\{c:c_i=\max_jc_j\},
 \qquad u=e_i\quad\hbox{on }\mathcal C_i.              \tag{7.7}
\]

The phase matrix is \(I_m/m\).  A continuous centroid field of magnitude
\(\gtrsim t^{-1/2}\) cannot jump from \(e_i\) to \(e_j\) across the tie
hyperplane: Lemma 1.1 forces a transition strip of physical width
\(\gtrsim\sqrt t\).

For independent standard Gaussian coordinates, the gap between the largest
and second-largest coordinates has scale \(1/\sqrt{\log m}\).  The
exact order-statistic formula is

\[
 \mathbb P(M_1-M_2>s)
 =m\int_{\mathbb R}\varphi(x)\Phi(x-s)^{m-1}dx.          \tag{7.8a}
\]

Mills' two-sided bounds, restricted to
\(\Phi(x)^{m}\in[1/4,3/4]\), give numerical \(c,C>0\) such that, for
\(m\ge3\) and \(0\le s\le c/\sqrt{\log m}\),

\[
 c s\sqrt{\log m}\le
 \mathbb P(M_1-M_2\le s)
 \le C s\sqrt{\log m}.                                 \tag{7.8b}
\]

For larger \(s\), the probability is bounded below by a numerical constant
and rapidly approaches one.  Taking
\(s\asymp\sqrt{t/P}=1/\sqrt{1+\tau}\) therefore gives

\[
 \nu\{\max c_i-\max_{j\ne i_*}c_j\lesssim\sqrt t\}
 \asymp
 \min\left\{1,{\sqrt{\log m}\over\sqrt{1+\tau}}\right\}.
                                                               \tag{7.8}
\]

Thus at \(m\asymp K/\tau\) and bounded \(\tau\), a numerical fraction of
the tilt law lies in the switching/tie region.  This is precisely the
separated-cell situation detected by Theorem 3.1; unlike the radial model,
there is no diffuse angular continuum.

However, to declare the tie strip profile-bad using (2.5) still requires
the defect scale (4.6), and finding a \(\tau\) with sufficiently small
relative profile charge still requires (6.4).  Equation (6.6) shows that
additive scale selection does not provide it.  Product cells therefore
validate the capacity lemma but do not rescue the multiscale closure.

The maximum-cell model is a phase-geometry stress test, not an asserted
KLS counterexample.  Product exponential, cube, simplex, and polyhedral
facet posteriors exhibit the same tie-region geometry; their locally sharp
phases arise only where the quadratic posterior core is separated from the
hard or one-sided boundary.

## 8. Formal conclusion

The following statements are now rigorous and dimension free:

1. the natural tilt law has the Poincare and Cheeger scales (1.3)--(1.4);
2. central separated phases have physical tilt separation (2.8);
3. fixed-mass separated phase cells force the capacity and profile charges
   (3.2)--(3.5);
4. a globally good phase law is coherent unless
   \((1+\tau)\Omega(e)\gtrsim1\);
5. the audited modulus leads to the explicit threshold (4.6).

They do not imply KLS.  High Fisher effective rank may be diffuse rather
than cell-like, conditioning on the good phase set has no available
dimension-free expansion theorem, noncentral tail corridors evade the
profile weight, and every polynomial defect threshold is integrable under
the exact logarithmic profile budget.  The mandatory examples show that
each loss is structural.

A viable continuation must prove both of the following new statements:

* a quantitative inverse theorem turning every diffuse low-defect phase law
  into a single common Gaussian factor of the original posterior family;
* a non-logarithmic profile or endpoint charge which can select the required
  factorization scale.

Neither statement follows from the pointwise angular estimate or from
ordinary Poincare expansion of the natural tilt law.  Consequently additive
\(\tau\)-scales do not beat the known pigeonhole mismatch.
