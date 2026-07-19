# Binary Fisher information versus the original calibrated rays

## 1. Scope and verdict

Let a balanced Kantorovich cut be disintegrated over its **original**
calibrated rays, and then expose the measure to the ordinary Gaussian
localization channel.  There is an exact posterior decomposition which
separates the desired normal alignment from one, and only one, obstruction:
posterior information about the label carried by the ray identity.

More precisely, at a posterior state the sharp-centroid defect controls the
weighted angular distance between the active Fisher direction and the
original ray normals, up to

\[
        \zeta=\operatorname {Var}\{\mathbb P(B=1\mid C,Q)\mid C\}.
\]

The quantity \(\zeta\) is not an informal error.  It is comparable to the
conditional mutual information \(I(B;Q\mid C=c)\), and its average is the
exact gap

\[
 I(B;C\mid Q)-I(B;C).
\]

Thus the bridge has a rigorous dichotomy:

1. small ray-identity information gives posterior-weighted normal
   alignment and long calibrated rays;
2. large ray-identity information is precisely the adaptive
   endpoint/rematching branch.

The Fisher operator bound gives pairwise near-orthogonality of the active
directions.  It transfers to the ray normals in branch 1.  What is not
proved is that branch 2 has positive-density cross-calibration.  A posterior
may select a different, very small ray packet for each observation.  The
finite Bregman inequality detects this only after paying the inverse
likelihood overlap of that packet.  No dimension-free bound for that overlap
is available.

## 2. Exact posterior reweighting of the ray quotient

Let

\[
 d\mu(x)=\int d\nu_q(r)\,d\eta(q),\qquad
 x=z_q+rN_q,\qquad |N_q|=1,                         \tag{2.1}
\]

be the nonbranching calibrated-ray disintegration of a signed-distance
Kantorovich potential \(f\), with \(f(z_q+rN_q)=r\).  Put

\[
 B={\bf1}_{\{r>0\}},\qquad
 \nu_q(r>0)=p\in(0,1)                                \tag{2.2}
\]

for \(\eta\)-almost every ray.  In particular \(B\) and \(Q\) are
independent under the original law.

For \(t>0\), observe

\[
 C=tX+\sqrt t\,G.                                    \tag{2.3}
\]

Conditionally on \(C=c\), the posterior is

\[
 d\pi_c(x)={1\over Z(c)}
  \exp\{c\cdot x-t|x|^2/2\}\,d\mu(x).               \tag{2.4}
\]

The posterior quotient law is \(r_c(dq)=L_c(q)\eta(dq)\), where

\[
\begin{aligned}
 Z_q(c)&=\int \exp\{c\cdot(z_q+rN_q)
                  -t|z_q+rN_q|^2/2\}\,d\nu_q(r),\\
 L_c(q)&={Z_q(c)\over\int Z_{q'}(c)d\eta(q')} .       \tag{2.5}
\end{aligned}
\]

Conditionally on \((C=c,Q=q)\), the ray law is

\[
 d\nu_{q,c}(r)={1\over Z_q(c)}
 \exp\{c\cdot(z_q+rN_q)-t|z_q+rN_q|^2/2\}
 d\nu_q(r).                                          \tag{2.6}
\]

It is one-dimensional \(t\)-strongly log-concave.  Define

\[
\begin{aligned}
 b_q&=\nu_{q,c}(r>0),&g&=\mathbb E_{r_c}b_q,\\
 m_q&=\mathbb E_{\nu_{q,c}}X,&m&=\mathbb E_{\pi_c}X,\\
 d_q&=\mathbb E[r\mid r>0,q,c]
       -\mathbb E[r\mid r<0,q,c],\\
 k_q&=\operatorname {Cov}_{\nu_{q,c}}(B,r)
       =b_q(1-b_q)d_q.                                \tag{2.7}
\end{aligned}
\]

All formulas below use this posterior quotient.  Replacing it silently by
the original law \(\eta\) is invalid.

## 3. The posterior normal-alignment dichotomy

Let

\[
 v=\operatorname {Cov}_{\pi_c}(B,X),\qquad
 u={v\over|v|},\qquad
 \zeta=\operatorname {Var}_{r_c}(b_q),                \tag{3.1}
\]

and assume \(v\ne0\).  Write

\[
 \delta_c=\mathcal I(g)-\sqrt t\,|v|\ge0.             \tag{3.2}
\]

Here \(\delta_c\) is the dimensionless sharp-centroid defect; it equals
\(\sqrt t\) times the defect denoted \(\Delta\) in the localization
notes.

**Proposition 3.1 (exact angular/polarization dichotomy).**  At every
posterior state,

\[
\boxed{
 \big\{\mathcal I(g)-\mathbb E_{r_c}\mathcal I(b_q)\big\}
 +\mathbb E_{r_c}\!\left[
      \mathcal I(b_q)-\sqrt t\,k_q\right]
 +\mathbb E_{r_c}\!\left[
      \sqrt t\,k_q(1-|u\cdot N_q|)\right]
 \le \delta_c+\sqrt\zeta .}                           \tag{3.3}
\]

Consequently

\[
\boxed{
 \mathbb E_{r_c}\!\left[
   \mathcal I(b_q)(1-|u\cdot N_q|)\right]
 \le\delta_c+\sqrt\zeta,}                             \tag{3.4}
\]

and

\[
\boxed{
 \mathbb E_{r_c}\!\left[
   \mathcal I(b_q)\min_{\epsilon\in\{-1,1\}}
       |u-\epsilon N_q|^2\right]
 \le2(\delta_c+\sqrt\zeta).}                         \tag{3.5}
\]

**Proof.**  Conditional covariance gives the exact decomposition

\[
\boxed{
 v=\mathbb E_{r_c}[k_qN_q]
       +\operatorname {Cov}_{r_c}(b_q,m_q).}           \tag{3.5a}
\]

The law of total covariance and Brascamp--Lieb for \(\pi_c\) give

\[
 \operatorname {Cov}_{r_c}(m_q)\preceq
 \operatorname {Cov}_{\pi_c}(X)\preceq t^{-1}I.
\]

Therefore

\[
 \sqrt t\,\left|u\cdot
       \operatorname {Cov}_{r_c}(b_q,m_q)\right|
 \le\sqrt\zeta.                                      \tag{3.6}
\]

The sharp one-dimensional centroid inequality applied to (2.6) says

\[
                  \sqrt t\,k_q\le\mathcal I(b_q).    \tag{3.7}
\]

Taking the scalar product of (3.5a) with \(u\), then replacing each signed
projection by its absolute value, and using (3.6), gives

\[
 \sqrt t|v|\le
 \mathbb E\{\sqrt t\,k_q|u\cdot N_q|\}+\sqrt\zeta.
\]

Subtract this from \(\mathcal I(g)\).  The first brace in (3.3) is
nonnegative by concavity of the Gaussian profile, and the other two terms
are nonnegative by (3.7).  This proves (3.3).  Moreover,

\[
\begin{aligned}
 \mathcal I(b_q)(1-|u\cdot N_q|)
 \le{}&\{\mathcal I(b_q)-\sqrt t\,k_q\}\\
 &+\sqrt t\,k_q(1-|u\cdot N_q|),
\end{aligned}
\]

which proves (3.4).  Finally use

\[
 \min_{\epsilon\in\{-1,1\}}|u-\epsilon N_q|^2
 =2(1-|u\cdot N_q|).
\]

The conclusion is deliberately for the unoriented normal line.  Recovering
the positive orientation requires a separate bound on the rays with
\(u\cdot N_q<0\); it cannot be obtained by multiplying (3.7) by a negative
number.  \(\square\)

There is also an unweighted, quantitative posterior statement.  Fix
\(0<\gamma<1/2\) and suppose

\[
                  \gamma\le g\le1-\gamma.             \tag{3.8}
\]

Put \(i_\gamma=\min_{a\in[\gamma/2,1-\gamma/2]}
\mathcal I(a)>0\).  Chebyshev and (3.3)--(3.5) give, for every \(a>0\),

\[
\boxed{
 r_c\{|b_q-g|>\gamma/2\}
 \le {4\zeta\over\gamma^2},                          \tag{3.9}
}
\]

and

\[
\boxed{
 r_c\left\{|b_q-g|\le\gamma/2,\ 
       \min_{\epsilon\in\{-1,1\}}|u-\epsilon N_q|\ge a\right\}
 \le {2(\delta_c+\sqrt\zeta)\over i_\gamma a^2}.    \tag{3.10}
}

Thus small \(\delta_c+\sqrt\zeta\) gives a fixed posterior mass of
calibrated normal **lines** in an arbitrarily prescribed fixed projective
cap about \(u\).  On those rays, both signs have fixed posterior mass.

There is also an explicit calibrated-length extraction.  Apart from an
\(r_c\)-set of mass at most

\[
 {4\zeta\over\gamma^2}
 +{2(\delta_c+\sqrt\zeta)\over i_\gamma},              \tag{3.11}
\]

one has both

\[
 b_q\in[\gamma/2,1-\gamma/2],\qquad
 \sqrt t\,k_q\ge {i_\gamma\over2}.                     \tag{3.12}
\]

Since \(k_q=b_q(1-b_q)d_q\) and \(b_q(1-b_q)\le1/4\),

\[
\boxed{
 d_q\ge {2i_\gamma\over\sqrt t}.}                      \tag{3.13}
\]

The positive and negative conditional means lie in the two halves of the
same maximal calibrated ray.  Therefore its two-sided geometric support
has total length at least \(2i_\gamma/\sqrt t\).  At \(t=K^{-1}\), these
are genuine \(c_\gamma\sqrt K\)-segments, not merely large active
one-dimensional variance in an unrelated direction.

## 4. The obstruction is exactly conditional ray-identity information

At a fixed observation, conditional mutual information is

\[
 \mathsf S(c):=I(B;Q\mid C=c)
 =\mathbb E_{r_c}\!\left[
   b_q\log{b_q\over g}+(1-b_q)\log{1-b_q\over1-g}
 \right].                                             \tag{4.1}
\]

Pinsker's inequality and the elementary chi-square upper bound give

\[
\boxed{
             2\zeta\le\mathsf S(c)
             \le {\zeta\over g(1-g)}.}                \tag{4.2}
\]

The upper bound follows from
\(D(\operatorname {Bern}(b)\|\operatorname {Bern}(g))
 \le(b-g)^2/[g(1-g)]\).

Since (2.2) says \(I(B;Q)=0\), the chain rule yields the exact averaged
identity

\[
\boxed{
 \mathbb E_C\mathsf S(C)
 =I(B;C\mid Q)-I(B;C).}                               \tag{4.3}
\]

Equivalently, in squared-error form,

\[
\boxed{
 \mathbb E_C\zeta(C)
 =\mathbb E\operatorname {Var}(B\mid C)
  -\mathbb E\operatorname {Var}(B\mid C,Q).}          \tag{4.4}
\]

Thus \(\zeta\) measures the improvement in label prediction obtained by
revealing the original ray after the Gaussian observation.  It can consume
a fixed part of the one-bit budget.  Neither (4.3) nor the scalar
I--MMSE identity makes it small.

## 5. Constant information forces fixed mass of genuinely long rays

The long-ray conclusion does not follow merely from a large average second
moment: that would permit a vanishing mass at an enormous scale.  The
binary entropy cap removes this escape.

Conditioned on \(Q=q\), the components of \(C\) perpendicular to \(N_q\)
are independent of \(B\).  The sufficient scalar observation is

\[
       {\langle C-tz_q,N_q\rangle\over\sqrt t}
       =\sqrt t\,T+G_1.                               \tag{5.1}
\]

Let \(\sigma_q^2=\operatorname {Var}_{\nu_q}T\).  Gaussian capacity and
data processing give

\[
\begin{aligned}
 I(B;C)&\le I(B;C,Q)
      =\int I(B;\sqrt tT+G_1\mid Q=q)d\eta(q)\\
 &\le\int\min\left\{h(p),{t\sigma_q^2\over2}\right\}d\eta(q).
                                                               \tag{5.2}
\end{aligned}
\]

**Proposition 5.1 (fixed-mass long-ray extraction).**  If
\(I(B;C)\ge i_0>0\), then

\[
\boxed{
 \eta\left\{q:\sigma_q^2\ge{i_0\over t}\right\}
 \ge {i_0\over2h(p)}.}                                \tag{5.3}
\]

**Proof.**  On the complement of the displayed set, the integrand in
(5.2) is at most \(i_0/2\); on the set it is at most \(h(p)\).  Therefore
\(i_0\le i_0/2+h(p)\eta(F)\).  \(\square\)

At \(t=K^{-1}\), (5.3) gives a fixed quotient mass of rays of standard
deviation at least \(\sqrt{i_0K}\).  In particular the conclusion is not
an exponentially small tilt-tail statement.

A Fisher lower bound at one isolated time does not imply the hypothesis of
Proposition 5.1.  A multiplicative phase does.  With the notation of
`binary_fisher_constraints.md`, if

\[
                 \operatorname {tr}R_s\ge r_0
 \quad(t/2\le s\le t),                                \tag{5.4}
\]

then the exact binary I--MMSE identity gives

\[
 I(B;C_t)-I(B;C_{t/2})
 =\int_{t/2}^{t}{\operatorname {tr}R_s\over2s}ds
 \ge {r_0\log2\over2}.                               \tag{5.5}
\]

Thus (5.3) applies with
\(i_0=r_0\log2/2\).

There is an independent extraction directly from a near-extremal T3
witness.  Let \(A=\mathbb E|f-\mathbb Ef|\), translate so that
\(\mathbb Ef=0\), and suppose \(\operatorname {Var}f\le K\).  If

\[
 d_q=\mathbb E[T\mid T>0,q]-\mathbb E[T\mid T<0,q],
\]

then exact ray balance gives

\[
 \mathbb E_\eta d_q={A\over2p(1-p)},\qquad
 \mathbb E_\eta d_q^2\le {K\over p(1-p)}.             \tag{5.6}
\]

The second inequality uses
\(\sigma_q^2\ge p(1-p)d_q^2\) and total variance.
Paley--Zygmund therefore gives

\[
\boxed{
 \eta\left\{d_q\ge {A\over4p(1-p)}\right\}
 \ge {A^2\over16p(1-p)K}.}                            \tag{5.7}
\]

Every ray in this set has

\[
 \sigma_q\ge\sqrt{p(1-p)}d_q
             \ge {A\over4\sqrt{p(1-p)}}.              \tag{5.8}
\]

Hence \(A^2\asymp K\), with \(p\) bounded away from zero and one,
already gives fixed original quotient mass at scale \(\sqrt K\).

## 6. Fisher effective rank gives pairwise orthogonality

At a fixed time put

\[
 w(c)={\eta_c^2\mathcal I(g_c)^2\over g_c(1-g_c)},
 \qquad
 R_t=\mathbb E_C[w(C)u(C)u(C)^T],
 \qquad M=\operatorname {tr}R_t.                      \tag{6.1}
\]

The exact covariance inequality gives \(R_t\preceq tI\).  Normalize
\(d\lambda(c)=w(c)d\mathbb P_C(c)/M\).  For independent
\(C,C'\sim\lambda\),

\[
\boxed{
 \mathbb E_{\lambda\otimes\lambda}
       \langle u(C),u(C')\rangle^2
 ={\operatorname {tr}(R_t^2)\over M^2}
 \le {t\over M}.}                                    \tag{6.2}
\]

Thus if \(t=K^{-1}\) and \(M\ge m_0>0\), two Fisher-weighted states have
squared overlap at most \((m_0K)^{-1}\) on average.  If Proposition 3.1
supplies coupled normal lines with
\(\mathbb E\min_{\epsilon=\pm1}|N-\epsilon u|^2\le\varepsilon\), choose
\(\widetilde N=\operatorname {sgn}(u\cdot N)N\).  Then

\[
 \mathbb E\langle \widetilde N,\widetilde N'\rangle^2
 \le {3\over m_0K}+6\varepsilon.                     \tag{6.3}
\]

This is the precise pairwise near-orthogonality statement.  An effective
rank lower bound by itself is weaker than a deterministic collection of
\(K\) mutually orthogonal vectors; (6.2), not that informal assertion, is
what transfers to packet estimates.

There is no hidden loss when passing from posterior ray incidence back to
the original quotient.  The joint law is

\[
 d\Xi(c,q)=d\mathbb P_C(c)\,r_c(dq),
 \qquad \Xi_Q=\eta.                                    \tag{6.4}
\]

On central near-equality states the weight \(w(c)\) in (6.1) is bounded
above and below by positive constants depending only on the centrality and
near-equality parameters.  Hence a fixed \(\lambda\)-mass of states, each
carrying a fixed \(r_c\)-mass of the aligned long rays extracted in
(3.9)--(3.13), gives fixed \(\Xi\)-mass and therefore fixed original
\(\eta\)-mass.  What this produces is a positive-density **incidence law**
\((c,q)\) with dispersed \(u(c)\) and projectively close \(N_q\).  It does
not produce calibrated chords between rays attached to two different
states; that is the cross-calibration gap addressed in Section 7.

## 7. What extremality converts, and the likelihood-overlap loss

Suppose, for simplicity, that all original ray laws lie in one balanced
scale band \([s,A s]\), so the finite Bregman lower bound has the form

\[
 \int B_q(h(q)-\mathbb E_\eta h)d\eta(q)
 \ge {c\over s}\operatorname {Var}_\eta h             \tag{7.1}
\]

for \(\operatorname {osc}h\le c_1s\).  At a posterior state let
\(L_c=dr_c/d\eta\), and assume on the family under consideration that

\[
                         L_c(q)\le M.                  \tag{7.2}
\]

Choose

\[
 h_c(q)=\theta s\{b_q-\mathbb E_\eta b_q\},
 \qquad0<\theta\le c_1.                               \tag{7.3}
\]

Since

\[
 \zeta=\inf_a\int L_c(q)(b_q-a)^2d\eta(q)
 \le M\operatorname {Var}_\eta b_q,                  \tag{7.4}
\]

the exact finite height-defect inequality for a true T3 extremizer implies

\[
\boxed{
 2\int|r_{h_c}|d\mu
 \ge c\theta^2s\operatorname {Var}_\eta b_q
 \ge {c\theta^2s\over M}\zeta.}                      \tag{7.5}
\]

Here \(f-h_c(Q)+r_{h_c}\) is any globally one-Lipschitz lattice or
McShane realization of the ideal ray height.  Thus large synergy really
does force switching/endpoint defect, but only after division by the
posterior quotient likelihood bound \(M\).

There is no dimension-free bound for \(M\).  A high-dimensional Gaussian
observation can select one of exponentially many ray packets, each of
exponentially small original quotient mass.  Averaging over observations
returns \(\eta\), but it does not turn the joint law
\(\mathbb P(dc)r_c(dq)\) into the product law
\(\mathbb P(dc)\eta(dq)\) used by a deterministic quotient height.  This
is the same all-amplitude saturation recorded in
`finite_medial_competitor.md`: multilevel heights do not acquire a packet
entropy factor.

The boundary--quotient lemma in `ray_theta_averaging.md` gives a
complementary exact estimate.  In heat coordinates \(s=t^{-1}\), with
\(W=\nabla P_s(B\mu)/P_s\mu\), it proves

\[
 \mathbb E\!\left[|W|
 \|uu^T-N_QN_Q^T\|_{HS}^2\right]
 \le12\delta_J(s)+16\mathfrak B_s,                    \tag{7.6}
\]

where

\[
 \mathfrak B_s=\int q_q(0)
     \int|T_sB-B|d\nu_q\,d\eta(q).                   \tag{7.7}
\]

Equation (7.6) is a direct active-direction-to-ray-normal bridge.  Its
residual \(\mathfrak B_s\) is another exact form of the same localization
problem: it mixes boundary density with raywise prediction error and does
not localize to a chosen long-ray band.  If every ray has scale at least
\(L\), then \(\mathfrak B_s\le2U(s)/L\); a merely positive-mass long band
does not give this bound.

## 8. A log-concave calibrated model with genuine synergy

Raywise balance does not make \(\zeta\) identically zero.  There is an
explicit globally log-concave, isotropic, signed-distance calibrated model
with a central posterior state and a fixed positive synergy.

In \(\mathbb R^2\), let

\[
 d\mu(u,v)\propto
 \exp\left[-{1\over3}(|u|+|v|)^2\right]\,du\,dv,
 \qquad f(u,v)={|u|-|v|\over\sqrt2}.                  \tag{8.1}
\]

The isotropy calculation is the \(m=1\) case of the Clifford model.  In
each quadrant, write

\[
 (u,v)={1\over\sqrt2}ig(s_u(R+r),s_v(R-r)\big),
 \quad R>0,\quad |r|<R.                               \tag{8.2}
\]

The ray normal is \((s_u,-s_v)/\sqrt2\), and the original conditional ray
law is uniform on \((-R,R)\).  It is exactly balanced.  Reflection
\(r\mapsto-r\) gives a calibrated coupling of the two sign restrictions,
so \(f\) is a Kantorovich potential for this cut.

Take posterior curvature \(t=1\) and natural tilt \(c=(1,1)\).  Put

\[
 Z_\beta(R)=\int_{-R}^R e^{\beta r-r^2/2}dr,
 \qquad
 b_\beta(R)={\int_0^R e^{\beta r-r^2/2}dr\over Z_\beta(R)}. \tag{8.3}
\]

The four posterior quotient weights, up to one common normalization, are

\[
\begin{array}{c|c|c}
(s_u,s_v)&\beta&\text{weight in }dR\\ \hline
(+,+)&0&e^{-7R^2/6+\sqrt2R}Z_0(R)\\
(+,-)&\sqrt2&e^{-7R^2/6}Z_{\sqrt2}(R)\\
(-,+)&-\sqrt2&e^{-7R^2/6}Z_{-\sqrt2}(R)\\
(-,-)&0&e^{-7R^2/6-\sqrt2R}Z_0(R).
\end{array}                                            \tag{8.4}
\]

Block swap preserves the posterior and flips the label, so \(g=1/2\).
The two same-sign quadrant families have \(b=1/2\), while the two cross
families have probabilities \(b_{\sqrt2}(R)\) and
\(1-b_{\sqrt2}(R)\).  Hence \(\zeta>0\).

The active centroid is nonzero and lies on the antisymmetric line
\(\mathbb R(1,-1)\).  Indeed block swap sends \(v\) to \(-v\).  Its scalar
product with \((1,-1)\) is positive: the two cross quadrants contribute
positively, while the \(++\) and \(--\) contributions have opposite signs
and the former has the strictly larger weight
\(e^{\sqrt2R}\) versus \(e^{-\sqrt2R}\).  Thus this is a genuine central
Fisher state, not a state at which the active direction is undefined.

For a completely explicit lower bound, restrict the two cross quadrants to
\(1\le R\le11/10\).  On this interval,

\[
 b_{\sqrt2}(R)-{1\over2}
 \ge {e^{-1/2}\sinh(1/\sqrt2)
       \over4(11/10)\cosh(11\sqrt2/10)}>0.04.          \tag{8.5}
\]

Their unnormalized mass on the interval is at least
\(0.1e^{-847/600}\), while the total normalizer is at most

\[
 8e^{6/7}\int_0^\infty R e^{-7R^2/12}dR<17.           \tag{8.6}
\]

Consequently

\[
                         \boxed{\zeta>2\cdot10^{-6}.}  \tag{8.7}
\]

This is a universal positive number, not an exponentially small
tilt-overlap assertion.  The example is not a large-scale T3 extremizer;
as in the higher-dimensional Clifford model, its long-scale portion escapes
to a Gaussian radial tail.  It proves only the needed logical point:
log-concavity, exact ray balance, signed-distance calibration, and a central
posterior do not set the synergy term to zero.  Fixed-mass large-scale
exclusion must use global extremality or a new endpoint theorem.

## 9. The convex-nullity proposal: precise target and present failure

Suppose branch 1 of Proposition 3.1 has transferred the Fisher direction
law to a fixed positive-mass, single-scale family of original normals, and
(6.2) has made those normals pairwise dispersed.  Long calibrated rays give
the pointwise shape bound

\[
                         s^2\|S_y\|_{HS}^2\le C.       \tag{9.1}
\]

This does **not** currently permit a classification by convex nullity.
The zero separator of a signed-distance extremizer need not be convex, its
shape operator may have both signs, and a long normal ray does not make a
relative-nullity leaf complete.  The Clifford cone realizes all three
warnings.  The identities

\[
 \nabla_TS=SC_T,\qquad SC_T=C_T^*S,
 \qquad\nabla_TC_T=C_T^2                             \tag{9.2}
\]

force a cylinder from complete nullity leaves only when \(S|_{\ker(S)^\perp}\)
is definite.  Mixed signs allow nontrivial splitting tensors.

A sufficient new statement is the following formal target.

> **Fisher-weighted focal-coherence lemma.**  Fix numerical
> \(\alpha,\delta,A>0\).  Let a \(C^3\) signed-distance separator of an
> isotropic log-concave law have a ray family \(F\) of quotient mass at
> least \(\alpha\), with balance in \([\delta,1-\delta]\) and
> \(s\le\sigma_y\le As\).  Suppose a Fisher-weighted coupling assigns to
> each ray a unit vector \(u_y\) such that
> \(\mathbb E_F|u_y-N_y|^2\le\varepsilon\) and the coupled direction
> covariance has operator norm at most \(Cs^{-2}\).  If an active
> curvature block has rank \(k\), energy
> \(\operatorname {tr}S_y^2\asymp s^{-2}\), and
> \(\|S_y\|_{op}^2\lesssim s^{-2}/k\) on a positive-mass subfamily, then
> either:
> 1. the normals and nullity leaves lie, up to \(o_{\varepsilon}(1)\), in
>    one fixed \(O(k)\)-dimensional cylindrical projection; or
> 2. there are a fixed rank-\(O(k)\) projection \(P\), a point \(z\), and
>    a positive-mass subfamily of focal leaves with
>    \(|P(Z_y-z)|\ge c s\sqrt k\) and
>    \(|Pz|\le(c/4)s\sqrt k\).

Alternative 1 would be charged by the covariance/projection budget.
Alternative 2 is exactly the coherent radial feature excluded by translated
thin shell or projected Paouris (`active_endpoint_coercivity.md`, Lemma
11.3; `finite_medial_competitor.md`, Lemma 8.1).

The target is strictly stronger than (9.1), relative nullity, or the Fisher
effective-rank bound.  None of those inputs controls how focal centers from
different curvature leaves move.  Therefore the proposed
``cylinder versus concurrent/radial'' step is not presently a consequence
of convex nullity; it is the same fixed-mass focal-coherence lemma already
isolated by the endpoint analysis.

## 10. Exact remaining bridge

The proved chain is

\[
\begin{gathered}
 \text{Fisher phase on a multiplicative window}
 \Longrightarrow
 \text{fixed original mass of }\sqrt K\text{-rays},\\
 R_t\preceq K^{-1}I, \operatorname {tr}R_t\asymp1
 \Longrightarrow
 \mathbb E\langle u,u'\rangle^2\lesssim K^{-1},\\
 \delta_c+\sqrt{\zeta_c}\ll1
 \Longrightarrow
 \text{posterior-weighted }N_Q\approx u.
\end{gathered}                                         \tag{10.1}
\]

What is not proved is either

\[
 \mathbb E_{\text{Fisher}}\zeta_c=o(1),               \tag{10.2}
\]

or a dimension-free conversion of \(\zeta_c\gtrsim1\) into the robust
positive-density cross-calibration needed by Sections 8--9 of
`finite_medial_competitor.md`.  Equation (7.5) shows the exact loss in the
obvious extremality conversion: it is the posterior quotient likelihood
ratio \(M\), which can encode exponentially many adaptively selected ray
packets.  Removing that loss, or proving the Fisher-weighted focal-coherence
lemma in Section 9, is the load-bearing new statement.
