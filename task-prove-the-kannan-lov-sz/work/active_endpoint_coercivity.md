# Active-endpoint coercivity: a smooth Clifford-cone survivor

## 1. Verdict

The proposed active-endpoint conclusion is false under the literal
hypotheses in the question.  There is an explicit globally log-concave,
isotropic measure and a globally defined signed-distance function for which:

1. every localization ray is exactly balanced;
2. on a positive-boundary-weight tail family \(F_s\), every conditional
   standard deviation is at least \(s\);
3. the exact curvature bound is
   \[
     \sigma_y^2\|S_y\|_{HS}^2<1;
   \]
4. both focal endpoints lie
   \[
     \sqrt{2m+1}
   \]
   conditional standard deviations from the zero surface, and their
   conditional endpoint densities are exactly zero;
5. the actual smooth second-variation form, with its genuine normal-metric
   coefficients, already satisfies
   \[
     \mathcal J_F(h)\ge P_F\operatorname {Var}_{\nu_F}h
   \]
   for every smooth ray height \(h\); and
6. all \(s^5\) effective-rank, projection, and localized projection budgets
   hold.

No ideal heat bath or comparison with endpoint probabilities is used.  The
coercivity is supplied by the smooth curvature form itself.

The example has an important and sharp limitation.  Its long-ray family is a
far radial tail and
\[
  P_F\le C\exp(-c_ms^2),
\]
in fact with an exponent of order \(ms^2\).  Thus it refutes
“\(P_F>0\) implies active endpoints,” but it does not refute a theorem with
a dimension-free lower bound on the boundary fraction, or with positive
bulk mass in one variance band.  That strengthened statement remains open.

## 2. The measures that must not be confused

For a ray quotient \((\Omega,\eta)\), write
\[
 b_y=q_y(0),\qquad
 P_F=\int_F b_y\,d\eta(y),\qquad
 d\nu_F(y)=\frac{b_y}{P_F}\mathbf1_F(y)\,d\eta(y).       \tag{2.1}
\]
The exact centered boundary term from signed-distance variation is
\[
 \int_F b_y(h-\mathbb E_{\eta_F}h)^2\,d\eta
 =P_F\left(
   \operatorname {Var}_{\nu_F}h+
   (\mathbb E_{\nu_F}h-\mathbb E_{\eta_F}h)^2
 \right),                                               \tag{2.2}
\]
when both laws are normalized on \(F\).  In particular,
\(P_F\operatorname {Var}_{\nu_F}h\) is the weaker, correctly normalized
quantity appearing in the proposed endpoint spectral statement.

The counterexample below proves coercivity for this weaker quantity on every
radial tail.  On the full quotient it also controls the additional centering
term in (2.2), so it is an honest second-variation-stable signed-distance
model rather than merely a normalized-form construction.

## 3. An isotropic globally log-concave Clifford model

Fix \(m\ge4\), write
\[
 \mathbb R^{2m}=\mathbb R_u^m\oplus\mathbb R_v^m,
\]
and put
\[
 \kappa_m=\frac{m+1}{2(2m+1)}.
\]
Let \(\mu_m\) have density
\[
 d\mu_m(u,v)=Z_m^{-1}
 \exp\!\left[-\kappa_m(|u|+|v|)^2\right]\,du\,dv.       \tag{3.1}
\]
The function \((u,v)\mapsto(|u|+|v|)^2\) is convex, so \(\mu_m\) is
globally log-concave.

It is also isotropic.  To check the normalization exactly, set
\[
 S=|u|+|v|,\qquad A=\frac{|u|}{|u|+|v|}.
\]
Polar coordinates in the two blocks show that \(S\) and \(A\) are
independent, with
\[
 \kappa_m S^2\sim\operatorname {Gamma}(m,1),
 \qquad A\sim\operatorname {Beta}(m,m).                 \tag{3.2}
\]
Consequently
\[
 \mathbb ES^2=\frac m{\kappa_m},\qquad
 \mathbb EA^2=\frac{m+1}{2(2m+1)}.
\]
Rotational symmetry within the \(u\)-block gives
\[
 \operatorname {Var}(u_i)
 =\frac1m\mathbb E(S^2A^2)
 =\frac{m+1}{2\kappa_m(2m+1)}=1,                       \tag{3.3}
\]
and the same holds in the \(v\)-block.  All cross covariances vanish.

Define
\[
 f(u,v)=\frac{|u|-|v|}{\sqrt2}.                         \tag{3.4}
\]
This is the signed Euclidean distance to the Clifford cone
\[
 \Sigma=\{(u,v):|u|=|v|\}.                              \tag{3.5}
\]
Indeed, after fixing the angular directions, the distance in the
\((|u|,|v|)\)-quadrant to the diagonal is
\(\lvert |u|-|v|\rvert/\sqrt2\); changing the angular directions can only
increase the distance.  Thus \(f\) is globally one-Lipschitz and
\(|\nabla f|=1\) away from the cone and the two coordinate axes.

The block-swap symmetry sends \(f\) to \(-f\), so the global cut is balanced.
More importantly, the balance holds on every individual normal ray.

## 4. Exact rays, endpoint distances, and curvature

Parametrize the regular cone by
\[
 y=y(R,\xi,\zeta)
   =\frac R{\sqrt2}(\xi,\zeta),
 \qquad R>0,\quad \xi,\zeta\in S^{m-1}.
                                                               \tag{4.1}
\]
The oriented unit normal is
\[
 N_y=\frac1{\sqrt2}(\xi,-\zeta).                       \tag{4.2}
\]
The maximal calibrated normal ray is
\[
 x(R,\xi,\zeta,t)
 =y+tN_y
 =\frac1{\sqrt2}\big((R+t)\xi,(R-t)\zeta\big),
 \qquad -R<t<R,                                        \tag{4.3}
\]
and \(f(x(R,\xi,\zeta,t))=t\).

The two completed focal endpoints are
\[
 B=(\sqrt2R\,\xi,0),\qquad
 C=(0,\sqrt2R\,\zeta),                                  \tag{4.4}
\]
so both half-lengths equal \(R\), and
\[
 B-C=2RN_y.                                             \tag{4.5}
\]

The shape operator of \(\Sigma\) has eigenvalues
\[
 0,\quad
 \underbrace{R^{-1},\ldots,R^{-1}}_{m-1},\quad
 \underbrace{-R^{-1},\ldots,-R^{-1}}_{m-1}.             \tag{4.6}
\]
Hence
\[
 J_R(t)=\det(I+tS_y)
       =\left(1-\frac{t^2}{R^2}\right)^{m-1}.            \tag{4.7}
\]

The ambient potential is constant along every ray:
\[
 \kappa_m\big(|u(t)|+|v(t)|\big)^2
 =2\kappa_mR^2.                                         \tag{4.8}
\]
Therefore the exact conditional density is
\[
 q_R(t)
 =\frac{\left(1-t^2/R^2\right)^{m-1}}
        {R\,B(1/2,m)}\,\mathbf1_{(-R,R)}(t).             \tag{4.9}
\]
It is even, so each ray assigns mass \(1/2\) to each sign.  Put
\[
 d_m=\frac1{B(1/2,m)}.
\]
Then
\[
 b_R=q_R(0)=\frac{d_m}{R},\qquad
 \sigma_R^2=\operatorname {Var}_{q_R}T
 =\frac{R^2}{2m+1}.                                     \tag{4.10}
\]
At the focal endpoints,
\[
 \frac{q_R(R-)}{q_R(0)}
 =\frac{q_R(-R+)}{q_R(0)}=0.                            \tag{4.11}
\]
Thus the endpoint-decay hypothesis is satisfied in its strongest possible
form.

Finally, (4.6) and (4.10) give the exact polynomial-log-concave curvature
budget
\[
 \boxed{\quad
 \sigma_R^2\|S_y\|_{HS}^2
 =\frac{2(m-1)}{2m+1}<1.
 \quad}                                                  \tag{4.12}
\]

## 5. Exact quotient and boundary laws

Put
\[
 c_m=2\kappa_m=\frac{m+1}{2m+1}.                       \tag{5.1}
\]
The surface element in the coordinates (4.1) is, up to the constant angular
normalization,
\[
 d\mathcal H^{2m-1}_\Sigma
 =\left(\frac R{\sqrt2}\right)^{2m-2}
   dR\,d\omega(\xi)\,d\omega(\zeta).                    \tag{5.2}
\]
Since the ray normalizer is
\[
 Z_R=e^{-c_mR^2}R B(1/2,m)
\]
up to the common ambient normalizing constant, the ray-quotient probability
\(\eta\) factors into uniform angular laws and the radial density
\[
 d\eta_R(R)\propto
 R^{2m-1}e^{-c_mR^2}\,dR.                               \tag{5.3}
\]
Multiplication by \(b_R=d_m/R\) gives the boundary law.  If
\[
 Y=c_mR^2,\qquad \alpha=m-\frac12,
\]
then
\[
 Y\sim\operatorname {Gamma}(m,1)\quad\text{under }\eta,
 \qquad
 Y\sim\operatorname {Gamma}(\alpha,1)\quad\text{under }\nu.
                                                               \tag{5.4}
\]
The two sphere variables remain independent and uniform under both laws.

Let
\[
 P=\int b_R\,d\eta
\]
be the total boundary weight.  It is bounded above and below by numerical
constants.  One exact expression is
\[
 P=\sqrt{\frac2\pi}\,\sqrt{\kappa_m}\,
 \frac{\Gamma(m+1/2)\Gamma(m-1/2)}{\Gamma(m)^2},         \tag{5.5}
\]
and the assertion follows from the standard gamma-ratio bounds.

For \(s\ge2\), define the long-ray tail family
\[
 F_s=\{R\ge R_0\},\qquad
 R_0=s\sqrt{2m+1}.                                      \tag{5.6}
\]
Then \(\sigma_R\ge s\) on \(F_s\), and
\[
 P_{F_s}
 =P\,\mathbb P\!\left\{
   G_\alpha\ge (m+1)s^2
 \right\},\qquad G_\alpha\sim\operatorname {Gamma}(\alpha,1).
                                                               \tag{5.7}
\]
This quantity is strictly positive.  The gamma Chernoff bound gives
\[
 0<P_{F_s}
 \le C\exp\!\left[
  -\alpha\big(s^2-1-\log s^2\big)
 \right]
 \le C e^{-c m s^2}.                                    \tag{5.8}
\]
Equation (5.8) is the weight escape which prevents this example from
settling a strengthened positive-boundary-fraction theorem.

## 6. The actual smooth second-variation coefficients

We now compute the geometric form itself.  Since \(p=q=1/2\), the tail
kernel from the signed-distance second variation is
\[
 \beta_R(t)=
 \begin{cases}
  \displaystyle\frac12\int_t^Rq_R(u)\,du,&t>0,\\[6pt]
  \displaystyle\frac12\int_{-R}^tq_R(u)\,du,&t<0.
 \end{cases}                                             \tag{6.1}
\]
It is even.

For a ray height \(h=h(R,\xi,\zeta)\), the cone metric and (4.6) give
\[
\begin{aligned}
 \left|(I+tS)^{-1}\nabla_\Sigma h\right|^2
 ={}&|\partial_Rh|^2\\
 &+\frac2{R^2(1+t/R)^2}|\nabla_\xi h|^2
 +\frac2{R^2(1-t/R)^2}|\nabla_\zeta h|^2.              \tag{6.2}
\end{aligned}
\]
Define the three scalar coefficients obtained by integrating (6.2) against
\(\beta_R(t)\,dt\).  The radial one is
\[
\begin{aligned}
 A_0(R)
 &=\int_{-R}^R\beta_R(t)\,dt
   =\frac12\mathbb E_{q_R}|T|\\
 &=\frac{Rd_m}{2m}.                                     \tag{6.3}
\end{aligned}
\]

The two angular coefficients agree by symmetry.  They can be evaluated
without any focal-model guess.  The exact one-dimensional trace identity is
\[
 \int_{-R}^R\beta_R(t)
 \frac{m-1}{R^2}
 \left[
  \frac1{(1+t/R)^2}+\frac1{(1-t/R)^2}
 \right]dt=b_R.                                         \tag{6.4}
\]
There is no ambient-convexity deficit because (4.8) is constant in \(t\),
and both endpoint terms vanish by (4.11).  Comparing (6.4) with (6.2)
therefore gives
\[
 A_\xi(R)=A_\zeta(R)=\frac{b_R}{m-1}.                   \tag{6.5}
\]

Consequently, for any radial family \(F\), the actual smooth form
\(\mathcal S_F\) satisfies the exact normalized identity
\[
\boxed{
 \frac{\mathcal S_F(h)}{P_F}
 =\mathbb E_{\nu_F}\!\left[
  \frac{R^2}{2m}|\partial_Rh|^2
  +\frac1{m-1}\left(
    |\nabla_\xi h|^2+|\nabla_\zeta h|^2
  \right)
 \right].}                                              \tag{6.6}
\]

This is the real normal-graph form.  It is not a heat bath and it contains
no endpoint conditional-variance replacement.

The normal chart becomes focal at \(t=\pm R\).  Formula (6.6) is nevertheless
the limit of the classical form on
\(|t|\le(1-\varepsilon)R\).  Near a focal endpoint,
\(\beta_R(t)\) vanishes to order \(m\), while the worst inverse metric in
(6.2) has order \((R-|t|)^{-2}\); hence the integrand is locally of order
\((R-|t|)^{m-2}\).  It is integrable for \(m\ge2\).  Any additional focal
epi-derivative is a nonnegative cost on the right side of stability.
Therefore the full smooth-plus-focal form obeys
\[
 \mathcal J_F(h)\ge\mathcal S_F(h).                     \tag{6.7}
\]
No comparison of focal conductance with endpoint probability is being made.

## 7. Exact stability of the smooth form

### 7.1 Boundary-variance stability on every radial tail

The boundary law on \(F_s\) is a product of the two uniform sphere laws and
the gamma law in (5.4), conditioned on
\[
 Y\ge y_0:=(m+1)s^2.                                    \tag{7.1}
\]
The product-sphere Poincare inequality is
\[
 \operatorname {Var}(h\mid R)
 \le\frac1{m-1}\mathbb E\!\left[
  |\nabla_\xi h|^2+|\nabla_\zeta h|^2\mid R
 \right].                                               \tag{7.2}
\]

For the radial part, write \(g(Y)=\mathbb E[h\mid R]\).  The gamma potential
on any interval in \((0,\infty)\) is
\[
 U(Y)=Y-(\alpha-1)\log Y+\text{constant},
 \qquad U''(Y)=\frac{\alpha-1}{Y^2}.                    \tag{7.3}
\]
The one-dimensional Brascamp--Lieb inequality, also valid after restriction
to the convex interval \([y_0,\infty)\), gives
\[
 \operatorname {Var}_{\nu_{F_s}} g
 \le\frac1{\alpha-1}
   \mathbb E_{\nu_{F_s}}\!\left[Y^2(g'(Y))^2\right].
                                                               \tag{7.4}
\]
On the other hand, since \(Y=c_mR^2\),
\[
 \frac{R^2}{2m}|\partial_Rg|^2
 =\frac2mY^2|g'(Y)|^2.                                  \tag{7.5}
\]
For \(m\ge3\),
\[
 \frac2m\ge\frac1{\alpha-1}
 =\frac1{m-3/2}.                                        \tag{7.6}
\]

Decompose total variance into conditional angular variance and variance of
the angular mean, use (7.2)--(7.6), and apply Jensen to the radial derivative
of that mean.  Equation (6.6) yields
\[
 \mathcal S_{F_s}(h)
 \ge P_{F_s}\operatorname {Var}_{\nu_{F_s}}h.           \tag{7.7}
\]
Combining this with (6.7) proves the exact stability postulated in the
question:
\[
 \boxed{\quad
 \mathcal J_{F_s}(h)
 \ge P_{F_s}\operatorname {Var}_{\nu_{F_s}}h
 \quad\text{for every smooth ray height }h.
 \quad}                                                  \tag{7.8}
\]

### 7.2 The full quotient also satisfies the exact centering inequality

The preceding tail argument proves precisely the normalization requested in
the active-endpoint statement.  There is also a stronger global check.
Under the full boundary law, \(Y\sim\operatorname {Gamma}(\alpha,1)\);
under the full ray quotient, \(Y\sim\operatorname {Gamma}(m,1)\).  Thus
\[
 w(Y):=\frac{d\eta}{d\nu}(Y)
 =\frac{\Gamma(\alpha)}{\Gamma(m)}Y^{1/2},
 \qquad \mathbb E_\nu w=1.                              \tag{7.9}
\]
Wendel's gamma-ratio inequality gives
\[
\begin{aligned}
 \operatorname {Var}_\nu w
 &=\alpha\left(\frac{\Gamma(\alpha)}{\Gamma(\alpha+1/2)}
        \right)^2-1\\
 &\le\frac1{2\alpha}=\frac1{2m-1}.                      \tag{7.10}
\end{aligned}
\]
Hence, for the radial angular mean \(g\),
\[
 (\mathbb E_\nu g-\mathbb E_\eta g)^2
 =\operatorname {Cov}_\nu(g,w)^2
 \le\frac1{2m-1}\operatorname {Var}_\nu g.              \tag{7.11}
\]

Equations (7.4)--(7.6), now without truncation, give the sharper factor
\[
 \mathbb E_\nu\!\left[
  \frac{R^2}{2m}|g_R|^2
 \right]
 \ge\left(2-\frac3m\right)\operatorname {Var}_\nu g.
                                                               \tag{7.12}
\]
For \(m\ge4\),
\[
 2-\frac3m\ge1+\frac1{2m-1}.                            \tag{7.13}
\]
Combining (7.2) and (7.9)--(7.13) proves
\[
 \mathcal S(h)
 \ge\int b_R\big(h-\mathbb E_\eta h\big)^2\,d\eta.      \tag{7.14}
\]
Thus the model satisfies the original signed-distance second-variation
stability inequality, including the exact \(\eta\)-centering, using the
actual smooth form alone.  This does not assert that \(f\) is the global
first-moment maximizer; it proves the precise local stability condition
under discussion.

## 8. The \(s^5\) budgets all hold

The normal law is independent of \(R\) under both \(\eta\) and \(\nu\), and
\[
 N=\frac1{\sqrt2}(\xi,-\zeta),\qquad
 \mathbb E_{\nu_{F_s}}NN^T=\frac1{2m}I_{2m}.            \tag{8.1}
\]
Thus the normal support rank and covariance effective rank are both \(2m\).
For every rank-\(d\) orthogonal projection \(Q\),
\[
 \mathbb E|QN|^2=\frac d{2m},\qquad
 \mathbb E|QN|^4\le\mathbb E|QN|^2=\frac d{2m}.         \tag{8.2}
\]
By (5.8), after adjusting numerical constants,
\[
 P_{F_s}s^3\le Cm,\qquad P_{F_s}s^5\le Cm.              \tag{8.3}
\]
Therefore
\[
 P_{F_s}s^3\mathbb E|QN|^2\le Cd,\qquad
 P_{F_s}s^5\mathbb E|QN|^4\le Cd,                       \tag{8.4}
\]
and
\[
 2m\ge cP_{F_s}s^5.                                     \tag{8.5}
\]
These are exactly the covariance, fourth-moment, and effective-rank budgets.

The localized projection budget also holds with the correct geometry.
If the projected endpoints \(QB,QC\) both lie in a ball of radius
\(\rho\), then (4.5) gives
\[
 |QN|\le\frac{\rho}{R}\le\frac{\rho}{R_0}.              \tag{8.6}
\]
When \(\rho\le R_0\), (8.3) and \(R_0^2=(2m+1)s^2\) imply
\[
\begin{aligned}
 P_{F_s}s^5\mathbb E\!\left[
 |QN|^4\mathbf1_{\{|QB-z|,|QC-z|\le\rho\}}
 \right]
 &\le P_{F_s}s^5\frac{\rho^4}{R_0^4}\\
 &\le C\rho^2.                                         \tag{8.7}
\end{aligned}
\]
For \(\rho\ge R_0\), use \(|QN|\le1\), (8.3), and
\(\rho^2\ge R_0^2\) to obtain the same bound.  Thus even the local
endpoint-ball version of the projection estimate sees no contradiction.

## 9. Failure of active-endpoint coercivity

Every ray in this model has
\[
 \frac{r_y^+}{\sigma_y}
 =\frac{r_y^-}{\sigma_y}
 =\sqrt{2m+1}.                                          \tag{9.1}
\]
Fix any proposed universal activity radius \(A\).  Choosing
\[
 m>\frac{A^2-1}{2}
\]
makes the set
\[
 \{y\in F_s:r_y^+\le A\sigma_y
              \text{ or }r_y^-\le A\sigma_y\}
\]
empty.  Its endpoint contribution is therefore zero for every height.
Nevertheless (7.8) gives full variance coercivity.  In fact (6.6) and
Section 7 show that the smooth form alone supplies it.

This disproves any assertion, based only on the hypotheses stated in the
question, that a universal positive fraction of the stability control must
come from endpoints at distance \(O(\sigma_y)\).  It also explains why a
spectral truncation in the normal coordinates cannot repair the argument:
the smooth curvature operator and the normal covariance are simultaneously
diagonal in this model.  On the two sphere factors,
\[
 \frac1{m-1}(-\Delta_{S^{m-1}})
\]
has spectral gap exactly one, while
\(\mathbb E NN^T=I/(2m)\).  Increasing the active normal rank merely adds
more first spherical harmonics on which the actual smooth form continues to
saturate boundary variance.

The local curvature estimate is likewise sharp in the relevant sense.
Although every principal curvature is only \(1/R\), there are
\(2m-2\) of them, and (4.12) remains of order one.  A dimension-free
Hilbert--Schmidt bound does not make the smooth form a small operator
perturbation.

## 10. What remains true under a nonnegligible-weight hypothesis

The survivor escapes through (5.8).  Consequently the following stronger
question is not answered by this construction:

> Suppose, in addition, that \(F\) carries a fixed fraction of boundary
> weight, or that it has fixed bulk mass in one variance band.  Must a fixed
> fraction of stability then come from \(O(\sigma)\)-endpoints?

The present inputs still do not prove that statement.  Two genuinely new
estimates would be needed:

1. an operator deficit for the smooth form, not merely
   \(\sigma^2\|S\|_{HS}^2\le C\); and
2. a comparison of the actual focal/medial coefficient
   \(\rho/|N_i-N_j|\) with the boundary-weighted endpoint incidence law.

The first requirement is exposed by the Clifford cone, where the smooth
form has gap one.  The second is exposed by the two-sheet coefficient
pathology: a small normal angle makes scalar focal conductance large without
making its translation trace large.  Neither comparison follows from the
\(s^5\) projection budgets.

Thus the exact status is:

* **proved:** the proposed active-endpoint coercivity is false when
  “positive boundary weight” means only \(P_F>0\);
* **proved:** the countermodel is globally log-concave, isotropic, globally
  signed-distance calibrated and balanced, satisfies the actual
  second-variation coefficients, and passes every stated projection budget;
* **not proved or refuted:** active-endpoint coercivity after imposing a
  dimension-free lower bound on boundary fraction or positive single-scale
  bulk mass.

That last, strengthened assertion is the remaining central-ray closure
problem.  It cannot be replaced by an ideal endpoint heat bath without the
missing coefficient comparison.

## 11. The sharpened fixed-mass claim: focal scale versus focal coherence

The counterexample identifies the right stronger question more precisely.
Assume now that a family lies in one scale band

\[
 s\leq \sigma_y\leq A_\sigma s,
 \qquad \eta(F)=\alpha>0.                              \tag{11.1}
\]

For a fixed balance parameter, the one-dimensional density estimates and
(2.1) give

\[
 P_F\asymp_{\delta,A_\sigma}\frac{\alpha}{s}.          \tag{11.2}
\]

Thus a fixed fraction under the boundary law and a fixed fraction under the
ray-quotient law are comparable on this band.  This is the regime in which
the exponentially small tail escape (5.8) is unavailable.

### 11.1 What an active curvature rank actually proves

The following elementary lemma is the precise pointwise content of the
suggested \(s\sqrt r\) focal scale.

**Lemma 11.1 (spectral focal-radius band).**  Let \(S\) be the shape
operator at the zero point of a balanced ray of scale \(\sigma\).  Suppose
that, for some integer \(k\geq1\),

\[
 c_0\sigma^{-2}\leq \Sigma:=\operatorname {tr}S^2
       \leq C_0\sigma^{-2},\qquad
 \operatorname {rank}S\leq Ak,qquad
 \|S\|_{op}^2\leq \frac{A\Sigma}{k}.                  \tag{11.3}
\]

Then there is a set \(I\) of at least \(k/(2A)\) principal curvatures such
that

\[
 \frac{c}{\sigma\sqrt k}
 \leq |\kappa_i|\leq
 \frac{C}{\sigma\sqrt k},qquad i\in I.               \tag{11.4}
\]

Equivalently, their formal focal points

\[
 z_i=y-\kappa_i^{-1}N_y                               \tag{11.5}
\]

lie at distances

\[
 c\sigma\sqrt k\leq |z_i-y|\leq C\sigma\sqrt k.      \tag{11.6}
\]

**Proof.**  Let

\[
 I=\left\{i:\kappa_i^2\geq
              \frac{\Sigma}{2Ak}\right\}.
\]

The complementary eigenvalues contribute at most \(\Sigma/2\), because
there are at most \(Ak\) of them.  The eigenvalues in \(I\) therefore carry
at least \(\Sigma/2\).  Since each is at most \(A\Sigma/k\),
\(|I|\geq k/(2A)\).  The two bounds on \(\kappa_i^2\), followed by the two
bounds on \(\Sigma\) in (11.3), give (11.4)--(11.6).  QED.

The relevant effective rank here is

\[
 r_{op}(S^2)=\frac{\operatorname {tr}S^2}{\|S\|_{op}^2},
                                                               \tag{11.7}
\]

together with comparability of support rank and effective rank.  The weaker
participation ratio

\[
 \frac{(\operatorname {tr}S^2)^2}{\operatorname {tr}S^4}
\]

only gives an a priori focal lower scale \(\sigma k^{1/4}\), not
\(\sigma\sqrt k\).  Likewise the Hilbert--Schmidt upper bound (4.12) alone
does not supply the lower bound in (11.3).  A rigorous use of “near
saturation” must therefore produce both a nontrivial curvature energy and
the operator-rank balance in (11.3); these cannot be silently folded into
the word *rank*.

In the Clifford model, \(k=2m-2\), all nonzero curvatures have equal
magnitude, and

\[
 R=\sigma_R\sqrt{2m+1}\asymp\sigma_R\sqrt k.           \tag{11.8}
\]

Thus the model attains the scale in Lemma 11.1.

### 11.2 The exact relative-nullity statement

There is a useful local cylinder theorem which does follow from Codazzi,
but it is weaker than the global cylinder assertion needed for closure.

**Lemma 11.2 (local relative-nullity leaves).**  On a \(C^3\) hypersurface
patch on which \(\operatorname {rank}S=k\) is constant, the distribution

\[
 \Delta=\ker S                                             \tag{11.9}
\]

is integrable.  Each of its leaves is an open subset of an affine
\((n-1-k)\)-plane, and the unit normal is constant on that leaf.

**Proof.**  If \(X,Y\in\Delta\), Codazzi gives

\[
 0=(\nabla_XS)Y-(\nabla_YS)X=-S[X,Y],                  \tag{11.10}
\]

so \(\Delta\) is integrable.  For any tangent \(Z\), self-adjointness and
Codazzi give

\[
\begin{aligned}
 \langle\nabla_XY,SZ\rangle
 &=-\langle Y,\nabla_X(SZ)\rangle\\
 &=-\langle Y,(\nabla_XS)Z\rangle
  =-\langle Y,(\nabla_ZS)X\rangle=0.
\end{aligned}                                           \tag{11.11}
\]

Thus \(\nabla_XY\in\Delta\).  Its ambient second fundamental part also
vanishes because \(\langle SX,Y\rangle=0\), so each leaf is ambient
totally geodesic.  Finally, \(D_XN=SX=0\).  QED.

If these affine leaves are parallel through a whole patch, the patch is a
generalized cylinder over a \(k\)-dimensional core, and all its normals lie
in one fixed \((k+1)\)-dimensional space.  This is the harmless low-core
model: taking the projection onto that space in the fourth-moment budget
gives, by (11.2),

\[
 \alpha s^4\leq Ck.                                    \tag{11.12}
\]

But Lemma 11.2 does **not** say that leaves belonging to different base
points are parallel.  The conullity splitting tensor can rotate the affine
leaves.  More precisely, with (E=\Delta^\perp), define

\[
 C_T X=-(\nabla_XT)_E,
 \qquad T\in\Delta,\quad X\in E.                       \tag{11.12a}
\]

For a parallel unit field (T) along a nullity geodesic, the Euclidean
Gauss--Codazzi equations give

\[
 \nabla_T S=SC_T,qquad SC_T=C_T^*S,qquad
 \nabla_TC_T=C_T^2.                                    \tag{11.12b}
\]

In a parallel frame the last equation has the explicit solution

\[
 C_T(t)=C_T(0)(I-tC_T(0))^{-1}.                         \tag{11.12c}
\]

Thus a two-sided nullity segment of length (2L) only rules out real
eigenvalues of (C_T(0)) larger than (1/L) in magnitude.  It does not
control nilpotent parts or nonreal spectrum.  If (S|_E) is definite,
(SC_T=C_T^*S) makes (C_T) self-adjoint in a positive definite metric;
complete two-sided leaves then force (C_T=0), recovering the cylindrical
case.  Mixed principal-curvature signs remove exactly this conclusion.

There is an even more basic issue here: a long *normal* ray does not make a
relative-nullity leaf complete; the two directions are perpendicular.  The
radial nullity leaves of the Clifford cone, which terminate at its vertex,
are an explicit warning against applying a global completeness/cylinder
theorem in the present hypotheses.

### 11.3 When Codazzi really does produce a center

There is an exact center statement for an umbilic block.  Suppose a
principal curvature \(\kappa\ne0\) has multiplicity at least two on a
constant-multiplicity patch, and let \(E_\kappa\) be its eigendistribution.
Codazzi implies

\[
 X\kappa=0\qquad(X\in E_\kappa).                       \tag{11.13}
\]

Indeed, choose a unit \(Y\in E_\kappa\) orthogonal to \(X\) and take the
inner product of
\((\nabla_XS)Y=(\nabla_YS)X\) with \(Y\).  The same equations make
\(E_\kappa\) integrable.  Along one of its leaves, define

\[
 z=y-\kappa^{-1}N.
\]

Then

\[
 D_Xz
 =X+\kappa^{-2}(X\kappa)N-\kappa^{-1}SX=0.          \tag{11.14}
\]

Thus the center is constant **on that curvature leaf**, which lies on a
sphere of radius \(|\kappa|^{-1}\).  Under (11.3), an umbilic active block
of dimension comparable to \(k\) consequently gives a \((k+1)\)-dimensional
radial feature of radius \(\asymp s\sqrt k\).

This contains the two canonical equality geometries.  A round generalized
cylinder has one such block.  The Clifford cone has two blocks, with
curvatures \(1/R\) and \(-1/R\); their focal maps are the two endpoint
spheres in the orthogonal \(u\)- and \(v\)-projections.

What (11.14) does not say is that centers of different curvature leaves
coincide, or even lie in one fixed low-dimensional projection.  Approximate
eigenvalue bands are weaker still: derivatives of their spectral
projections involve inverse eigenvalue gaps.  Neither scalar smooth
stability nor \(\operatorname {tr}S^2\) controls those derivatives.

### 11.4 Projected Paouris closes a coherent radial feature

For completeness, here is the exact probabilistic closure once the missing
geometric coherence has been obtained.

**Lemma 11.3 (coherent feature exclusion).**  Let \(Q\) have rank
\(d\asymp k\).  Suppose a subfamily \(G\) has quotient mass at least
\(\alpha_0\), all its conditional scales are comparable to \(s\), and for
some fixed \(z\in\operatorname {ran}Q\),

\[
 |Q(Z_y-z)|\geq a s\sqrt k\quad(y\in G),
 \qquad |Qz|\leq \frac a4s\sqrt k.                    \tag{11.15}
\]

Here \(Z_y\) is the zero point of the ray.  Then, for \(k\) larger than a
constant depending only on balance and \(a\),

\[
 \alpha_0\leq C\exp(-c_as\sqrt k).                     \tag{11.16}
\]

**Proof.**  Balance and one-dimensional log-concavity imply
\(|\mathbb E(T_y\mid y)|\leq C_\delta s\), and hence
\(\mathbb E(T_y^2\mid y)\leq C_\delta s^2\).  Chebyshev therefore gives,
uniformly on \(G\), a fixed positive conditional probability of

\[
 |T_y|\leq \frac a4s\sqrt k.
\]

On this event, \(|QN_y|\leq1\) and (11.15) imply

\[
 |QX|\geq c_as\sqrt k.                                 \tag{11.17}
\]

The marginal \(QX\) is isotropic and log-concave in dimension
\(d\asymp k\).  Paouris' deviation inequality applied to (11.17) gives
(11.16).  QED.

The center bound in (11.15) is not an artificial origin convention.  It
follows, for example, if the common center is within \(O(s)\) of the
conditional mean of \(QX\) on \(G\).  Jensen and isotropy then give

\[
 |Qz|\leq \sqrt{d/\alpha_0}+O(s),                      \tag{11.18}
\]

which is \(o(s\sqrt k)\) for fixed \(\alpha_0\) in the large-scale,
large-rank regime.  Thus one may derive center control from a genuinely
centered radial leaf; it need not be assumed to be the global barycenter.

Combining (11.12) and Lemma 11.3 closes a coherent generalized-cylinder
branch: a core too small violates the projection budget, while a large
umbilic core produces an \(s\sqrt k\) projected radial tail forbidden at
fixed mass.

### 11.5 The precise missing geometric step

The currently available hypotheses stop between Lemmas 11.1 and 11.3.
They give pointwise focal radii but not focal coherence:

* the smooth translation trace is invariant under rotation of the
  principal eigenspaces and records only scalar curvature charge;
* the \(s^5\) projection budgets control the ray normals \(N_y\), not the
  active tangent spaces or the focal-center maps;
* relative nullity gives one affine leaf at a time, not a common core
  projection; and
* Codazzi makes the center constant on an exact repeated-curvature leaf,
  not across distinct leaves or approximate spectral blocks.

Accordingly, the sharpened central-ray closure reduces to the following
concrete statement:

> From fixed single-scale mass, near-saturation of the actual smooth form,
> and a rank-balanced active curvature block, produce a fixed projection
> of rank \(O(k)\) and a positive-mass collection of curvature leaves whose
> centers are coherent in the sense of (11.15)--(11.18); or construct a
> globally log-concave signed-distance example in which these centers
> genuinely de-cohere.

The first alternative would finish by projected Paouris.  No fixed-mass
example realizing the second alternative is presently known.  The Clifford
model realizes the pointwise equality geometry, but isotropy pushes its
scale-\(s\) portion into the tail (5.8), exactly as Lemma 11.3 predicts.
Thus the literal active-endpoint claim is refuted, whereas the fixed-mass
claim is reduced to a specific Codazzi/global-coherence problem rather than
to an endpoint heat-bath comparison.
