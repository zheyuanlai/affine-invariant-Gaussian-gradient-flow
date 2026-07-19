# Binary--quotient synergy on calibrated rays

## Executive conclusion

There are two different sources of conditional dependence between the sign
label and the ray quotient after a Gaussian observation:

1. the geometric scalar tilt
   \(a_q(c)=\langle c-tz_q,N_q\rangle\) varies with the ray; and
2. even at one fixed scalar tilt, the balanced one-dimensional laws
   \(\nu_q\) may have different shapes.

The second source cannot be discarded.  Section 3 gives an explicit,
full-dimensional, isotropic log-concave countermodel in which all calibrated
directions are parallel and \(a_q(c)\) is independent of \(q\), but

\[
 \int \operatorname {Var}_{Q\mid C=c}
       \mathbb P(B=1\mid Q,C=c)\,d\mathbb P_C(c)>2\cdot10^{-21}.       \tag{0.1}
\]

For the common-center projection height

\[
 h_c^{\rm proj}(q)=\frac{a_q(c)}t,
\]

both the centered finite Bregman gain and the minimal switching defect are
exactly zero in this example.  Consequently no universal inequality can
bound binary--quotient synergy by the finite competitor associated only with
this projection height.  In fact the same is true for every height of the
form \(H_t(a_q(c))\), with \(H_t\) independent of the conditional ray law:
it is quotient-constant in the countermodel.  In particular, positive
synergy does not by itself imply angular dispersion, radial concurrence, or
Clifford geometry.

There is nevertheless an exact way to encode *all* synergy as a quotient
height.  If the balanced ray laws have density at zero at least \(\beta/L\),
let \(w_c=d\eta_c/d\eta\) be the posterior quotient density and put

\[
 h_c^{\rm syn}(q)=Lw_c(q)\{b_q(c)-g(c)\},
 \qquad g(c)=\int b_q(c)w_c(q)d\eta(q).                         \tag{0.2}
\]

This height is already centered under the original quotient law.  Its finite
Bregman gain \(\mathcal G_c\) satisfies

\[
 \mathcal G_c\ge {\kappa_\beta L\over4}\,\zeta(c)^2,
 \qquad
 \kappa_\beta=\min\{\beta/2,1/4\},                            \tag{0.3}
\]

where \(\zeta(c)=\operatorname {Var}_{\eta_c}b_q(c)\).  If the calibrated
potential is a true first-moment extremizer, its minimal finite switching
correction \(\mathcal D_f(h_c^{\rm syn})\) therefore obeys

\[
 \boxed{\quad
 \zeta(c)\le
 \left({8\mathcal D_f(h_c^{\rm syn})\over\kappa_\beta L}\right)^{1/2}.
 \quad}                                                       \tag{0.4}
\]

Thus order-one synergy on rays of scale \(L\) costs order \(L\) switching
mass.  The estimate uses no quotient expansion or overlap assumption.  It
does not by itself close the long-ray branch: the height in (0.2) contains
the possibly very concentrated likelihood ratio \(w_c\), and no geometric
upper bound on its switching correction follows from the Gaussian channel.
The countermodel proves why replacing (0.2) by the geometrically controlled
projection height loses essential information.

Nor can one average (0.2) over the observation.  Section 5 proves
\(\mathbb E_C h_C^{\rm syn}=L\,\mathbb E[(I-A^*A)B\mid Q]\); this is
identically zero in the diamond by label-reflection symmetry.  The only
general likelihood-free switching bound is of order \(L\), and a scaled
version of the same calibrated model shows that this order is sharp for the
channel identities alone.

## 1. Exact posterior formulas

Let

\[
 X=z_q+TN_q,
 \qquad d\mu(X)=d\nu_q(T)d\eta(q),
 \qquad B={\bf1}_{\{T>0\}},                                  \tag{1.1}
\]

where \(|N_q|=1\) and

\[
 \nu_q(T>0)=\nu_q(T<0)=\frac12                              \tag{1.2}
\]

for \(\eta\)-almost every \(q\).  Observe

\[
 C=tX+\sqrt tG,
 \qquad G\sim N(0,I),                                      \tag{1.3}
\]

independently.  The posterior likelihood, after deleting a factor depending
only on \(c\), is

\[
 \exp\{c\cdot x-t|x|^2/2\}.
\]

Consequently, conditional on \((Q,C)=(q,c)\), the ray law is

\[
 {e^{a_q(c)T-tT^2/2}d\nu_q(T)\over M_q(a_q(c))},
 \qquad
 a_q(c)=\langle c-tz_q,N_q\rangle,                         \tag{1.4}
\]

where

\[
 \begin{aligned}
 M_q(a)&=\int e^{aT-tT^2/2}d\nu_q(T),\\
 M_q^+(a)&=\int_{T>0}e^{aT-tT^2/2}d\nu_q(T),\\
 b_q(c)&={M_q^+(a_q(c))\over M_q(a_q(c))}.
 \end{aligned}                                             \tag{1.5}
\]

The posterior quotient density is

\[
 w_c(q):={d\eta_c\over d\eta}(q)
 ={e^{c\cdot z_q-t|z_q|^2/2}M_q(a_q(c))
   \over
   \int e^{c\cdot z_r-t|z_r|^2/2}M_r(a_r(c))d\eta(r)}.     \tag{1.6}
\]

In particular,

\[
 g(c)=\int w_c(q)b_q(c)d\eta(q),\qquad
 \zeta(c)=\int w_c(q)(b_q(c)-g(c))^2d\eta(q).              \tag{1.7}
\]

All these identities are direct disintegration identities; they require no
smoothness of the ambient density or of the ray quotient.

For reference, conditional mutual information dominates this variance.  In
natural logarithms, Pinsker's inequality for Bernoulli laws gives

\[
 I(B;Q\mid C=c)
 =\int w_c(q)\,D(\operatorname {Ber}(b_q(c))
                 \Vert\operatorname {Ber}(g(c)))d\eta(q)
 \ge2\zeta(c).                                           \tag{1.8}
\]

Thus a positive lower bound for the integrated \(\zeta\) is also a positive
lower bound for genuine \(B\)--\(Q\) synergy.

## 2. An exact synergy height and its finite switching cost

For a balanced ray law \(\nu_q\), define its centered finite Bregman
function

\[
 B_q(s)=\int|T-s|d\nu_q(T)-\int|T|d\nu_q(T).              \tag{2.1}
\]

Assume that every \(\nu_q\) has a log-concave density \(\varphi_q\) and,
on the ray family under consideration,

\[
                         \varphi_q(0)\ge {\beta\over L}.   \tag{2.2}
\]

Here \(L>0\) is the common ray scale and \(\beta>0\) is numerical.  The
one-dimensional finite Bregman estimate gives, for every real \(s\),

\[
 B_q(s)\ge {1\over2}
       \min\{\varphi_q(0)s^2,\tfrac12|s|\}
 \ge \kappa_\beta L
       \min\{(s/L)^2,|s|/L\},                            \tag{2.3}
\]

where

\[
                         \kappa_\beta=\min\{\beta/2,1/4\}.
\]

The following proposition is the quotient-weight normalization which is
missing from the naive projection height.

**Proposition 2.1 (posterior synergy as finite Bregman gain).**  Fix an
observation \(c\), abbreviate

\[
 d_q=b_q(c)-g(c),\qquad w_q=w_c(q),qquad
 h_q=Lw_qd_q.                                           \tag{2.4}
\]

Then \(h\in L^1(\eta)\), \(\int h_qd\eta(q)=0\), and

\[
 \boxed{\quad
 \mathcal G_c:=\int B_q(h_q)d\eta(q)
 \ge{\kappa_\beta L\over4}\,\zeta(c)^2.
 \quad}                                                  \tag{2.5}
\]

**Proof.**  Since \(|d_q|\le1\) and \(\int w_qd\eta=1\), the height is
integrable.  Equation (1.7) gives

\[
 \int h_qd\eta=L\int w_q(b_q-g)d\eta=0.                 \tag{2.6}
\]

Put \(y_q=w_q|d_q|\) and

\[
 A=\int\min\{y_q^2,y_q\}d\eta(q).
\]

By (2.3), \(\mathcal G_c\ge\kappa_\beta LA\).  On
\(\{y_q\le1\}\),

\[
 w_qd_q^2=y_q|d_q|\le y_q
 =\sqrt{\min\{y_q^2,y_q\}},                            \tag{2.7}
\]

whereas on \(\{y_q>1\}\),

\[
 w_qd_q^2=y_q|d_q|\le y_q
 =\min\{y_q^2,y_q\}.                                   \tag{2.8}
\]

Cauchy--Schwarz therefore yields

\[
 \zeta(c)=\int w_qd_q^2d\eta
 \le\sqrt A+A.                                         \tag{2.9}
\]

Since \(0\le\zeta(c)\le1/4\), (2.9) implies
\(A\ge\zeta(c)^2/4\): this is immediate if \(A\ge1\), and if
\(A<1\), then \(\sqrt A+A\le2\sqrt A\).  This proves (2.5). \(\square\)

For any calibrated potential \(f(X)=T\), and any integrable centered
quotient height, define

\[
 \mathcal D_f(h)=\inf\left\{
   \int|r|d\mu: f-h(Q)+r\text{ has a one-Lipschitz representative}
                         \right\}.                     \tag{2.10}
\]

The set in (2.10) is nonempty: \(r=h(Q)\) returns the original function
\(f\).  Suppose now, in addition, that \(f\) is a true maximizer of the
centered first-moment functional among one-Lipschitz functions.  The exact
finite height-defect inequality, applied along an
infimizing sequence, gives

\[
 \int B_q(h_q)d\eta(q)\le2\mathcal D_f(h).             \tag{2.11}
\]

Combining (2.5) and (2.11) proves

\[
 \boxed{
 \mathcal D_f(h_c^{\rm syn})
 \ge{\kappa_\beta L\over8}\zeta(c)^2,
 \qquad
 \zeta(c)
 \le\left({8\mathcal D_f(h_c^{\rm syn})
                 \over\kappa_\beta L}\right)^{1/2}.}  \tag{2.12}
\]

After integrating and applying Cauchy--Schwarz,

\[
 \boxed{\quad
 \int\zeta(c)d\mathbb P_C(c)
 \le\left\{ {8\over\kappa_\beta L}
       \int\mathcal D_f(h_c^{\rm syn})d\mathbb P_C(c)
       \right\}^{1/2}.
 \quad}                                                  \tag{2.13}
\]

For completeness, the integrand in (2.13) is measurable.  The space
\(L^1(\mu)\) is separable, so the subset of one-Lipschitz functions has a
countable \(L^1(\mu)\)-dense subfamily \((u_j)\).  Joint measurability of
\((c,q)\mapsto h_c^{\rm syn}(q)\), together with
\(\int|h_c^{\rm syn}|d\eta\le L\), makes
\(c\mapsto f-h_c^{\rm syn}(Q)\) strongly measurable as an
\(L^1(\mu)\)-valued map.  Therefore

\[
 \mathcal D_f(h_c^{\rm syn})
 =\inf_j\|u_j-(f-h_c^{\rm syn}(Q))\|_{L^1(\mu)}          \tag{2.13a}
\]

is measurable.  This also justifies taking the ordinary integral rather
than an outer integral in (2.13).

No absolute continuity bound on \(w_c\), quotient Poincare inequality, or
pairwise posterior overlap was used.  At a hypothetical long-ray scale
\(L\asymp\sqrt K\), (2.12) says that order-one synergy forces switching
cost of order \(\sqrt K\).  What is unavailable is an upper bound smaller
than this: \(h_c^{\rm syn}\) can inherit arbitrary quotient oscillation from
\(w_c\).  The trivial correction \(r=h_c^{\rm syn}(Q)\) only gives

\[
 \mathcal D_f(h_c^{\rm syn})
 \le L\int w_c|b_q-g|d\eta
 \le L\sqrt{\zeta(c)},                                 \tag{2.14}
\]

which is fully compatible with (2.12).

## 3. Isotropic log-concave countermodel to the projection-height mechanism

Let \(R=\sqrt6\) and let \(\mu\) be the uniform probability measure on the
two-dimensional diamond

\[
 K=\{(q,T)\in\mathbb R^2:|q|+|T|\le R\}.              \tag{3.1}
\]

The body is convex, so \(\mu\) is log-concave.  It is full-dimensional and
has barycenter zero.  Its area is \(2R^2=12\), and direct integration gives

\[
 \mathbb E q^2=\mathbb E T^2={R^2\over6}=1,
 \qquad \mathbb E(qT)=0.                               \tag{3.2}
\]

Thus \(\mu\) is isotropic.

Take

\[
 f(q,T)=T,qquad Q=q,qquad z_q=(q,0),qquad N_q=e_2.   \tag{3.3}
\]

The maximal calibrated ray through \(q\) is

\[
 T\in[-w(q),w(q)],\qquad w(q)=R-|q|,                   \tag{3.4}
\]

and the conditional law \(\nu_q\) is uniform on this interval.  Hence

\[
 B={\bf1}_{\{T>0\}}\quad\hbox{is independent of }Q.    \tag{3.5}
\]

Every ray direction is the same vector \(e_2\).

Use channel precision \(t=1\), so \(C=X+G\).  Given \(Q=q,C=c\), the
conditional density of \(T\) is proportional on \([-w(q),w(q)]\) to

\[
                         e^{c_2T-T^2/2}.                \tag{3.6}
\]

For \(a\in\mathbb R,w>0\), put

\[
 A_w(a)=\int_0^w e^{as-s^2/2}ds,qquad
 \beta_w(a)={A_w(a)\over A_w(a)+A_w(-a)}.              \tag{3.7}
\]

Then

\[
                         b_q(c)=\beta_{w(q)}(c_2).      \tag{3.8}
\]

For every \(a>0\), \(w\mapsto\beta_w(a)\) is strictly increasing.  Indeed,
writing \(Z=A_w(a)+A_w(-a)\), differentiation gives

\[
 {\partial\beta_w(a)\over\partial w}
 ={2e^{-w^2/2}\over Z^2}
   \int_0^w e^{-s^2/2}\sinh(a(w-s))ds>0.               \tag{3.9}
\]

For \(1\le a\le5/4\) and \(R/4\le w\le3R/4\), restricting the integral
in (3.9) to \([0,w/2]\), using \(\sinh(w/2)\ge w/2\), and using
\(Z\le2we^{aw}\), gives the completely explicit bound

\[
 {\partial\beta_w(a)\over\partial w}
 \ge {1\over8}
   \exp\left\{-{135\over64}-{15R\over8}\right\}.      \tag{3.10}
\]

Consequently

\[
 \beta_{3R/4}(a)-\beta_{R/4}(a)
 \ge d_0:={R\over16}
   \exp\left\{-{135\over64}-{15R\over8}\right\}.      \tag{3.11}
\]

Consider the observation rectangle

\[
 E=\{|c_1|\le1/4,\ 1\le c_2\le5/4\}.                 \tag{3.12}
\]

The two quotient bands

\[
 H_0=\{|q|\le R/4\},
 \qquad H_1=\{3R/4\le|q|\le7R/8\}                    \tag{3.13}
\]

have original probabilities

\[
                         \eta(H_0)={7\over16},qquad
                         \eta(H_1)={3\over64}.          \tag{3.14}
\]

For \(c\in E\), the log-likelihood
\(c\cdot x-|x|^2/2\) on \(K\) has oscillation at most

\[
 {61\over16}+{\sqrt{39}\over2}<7.                     \tag{3.15}
\]

It follows directly from the posterior density that

\[
 \eta_c(H_0)\ge {7\over16}e^{-7},
 \qquad
 \eta_c(H_1)\ge {3\over64}e^{-7}.                     \tag{3.16}
\]

On \(H_0\), \(w(q)\ge3R/4\), and on \(H_1\),
\(w(q)\le R/4\).  Equations (3.8)--(3.11), followed by the pairwise formula
for variance, therefore give, for every \(c\in E\),

\[
 \zeta(c)\ge {63\over131072}
 \exp\left\{-14-{135\over32}-{15R\over4}\right\}.     \tag{3.17}
\]

The event \(E\) itself has explicitly positive channel probability.  If

\[
 |X_1|,|X_2|\le1/16,\quad
 G_1\in[-3/16,3/16],\quad
 G_2\in[17/16,19/16],                                  \tag{3.18}
\]

then \(C\in E\).  The first event in (3.18) has \(\mu\)-probability
\(1/768\).  Lower-bounding the two Gaussian interval probabilities by
interval length times the minimum density gives

\[
 \mathbb P(C\in E)
 \ge {1\over32768\pi}\exp\{-185/256\}.                \tag{3.19}
\]

Multiplying (3.17) and (3.19) proves

\[
 \boxed{
 \int\zeta(c)d\mathbb P_C(c)
 \ge {63\over2^{32}\pi}
 \exp\left\{-{4849\over256}-{15\sqrt6\over4}\right\}
 >2\cdot10^{-21}.}                                     \tag{3.20}
\]

This tiny displayed number is chosen only to make every constant auditable;
the actual synergy is much larger.

The effect is not produced only by the short rays near the vertices.  It is
already witnessed on a positive-mass family whose conditional half-lengths
are all comparable.  Namely, put

\[
 J_0=\{|q|\le R/8\},\qquad
 J_1=\{R/4\le|q|\le3R/8\}.                             \tag{3.20a}
\]

Their original masses are \(15/64\) and \(11/64\), respectively, and every
ray in \(J_0\cup J_1\) has

\[
                         5R/8\le w(q)\le R.             \tag{3.20b}
\]

On \(J_0\), \(w\ge7R/8\), while on \(J_1\), \(w\le3R/4\).
Repeating (3.10) on \(3R/4\le w\le7R/8\) gives

\[
 \beta_{7R/8}(a)-\beta_{3R/4}(a)
 \ge {R\over64}
   \exp\left\{-{735\over256}-{35R\over16}\right\}.       \tag{3.20c}
\]

The same posterior oscillation estimate and the same event \(E\) therefore
give the independent explicit lower bound

\[
 \boxed{
 \int\zeta(c)d\mathbb P_C(c)
 \ge {495\over2^{38}\pi}
 \exp\left\{-{5239\over256}-{35\sqrt6\over8}\right\}
 >10^{-23},}                                            \tag{3.20d}
\]

using only the two comparable-scale ray bands in (3.20a).  Thus neither
vanishing ray mass nor a mixture of microscopic and macroscopic ray scales
causes the phenomenon.

On the other hand, the centered common-center projection height is exactly

\[
 h_c^{\rm proj}(q)
 ={\langle c-z_q,e_2\rangle}=c_2,                      \tag{3.21}
\]

independent of \(q\).  Hence

\[
 h_c^{\rm proj}-\eta h_c^{\rm proj}=0,
 \qquad
 \int B_q(h_c^{\rm proj}-\eta h_c^{\rm proj})d\eta=0. \tag{3.22}
\]

Thus the centered height to which (2.10) applies is zero and
\(\mathcal D_f(0)=0\).  Equivalently, before centering,
\(f-h_c^{\rm proj}(Q)=f-c_2\) is already one-Lipschitz.
Equations (3.20)--(3.22) disprove
any proposed estimate of synergy by the centered Bregman and switching
defects of \(a_q(c)/t\).  More generally, because \(a_q(c)=c_2\) for every
ray, any scalar-tilt-only height \(H_t(a_q(c))\) is quotient-constant and has
the same zero defects.

## 4. What the countermodel rules out, and what it does not

The diamond is a fixed-dimensional, fixed-scale model.  It does not produce
a bad KLS sequence and it does not refute a theorem which adds a quantitative
asymptotic homogeneity assumption on long conditional ray laws.  Although
\(f=T\) is a calibrated one-Lipschitz potential, no claim is made here that
it is a global \(T3\) extremizer for the diamond.  In fact it is not:
under the change of variables \(u=q+T,v=q-T\), the diamond becomes the
square \(\{|u|,|v|\le R\}\), and hence

\[
 \mathbb E\left|{q+T\over\sqrt2}\right|
 ={R\over2\sqrt2}>{R\over3}=\mathbb E|T|.              \tag{4.1}
\]

Thus the model refutes any
deduction made from calibrated disintegration, isotropy, log-concavity, and
the channel formulas alone; a new use of global extremality which forces
ray-law homogeneity could evade it.  It does
prove the following points without qualification:

1. Variation of \(a_q(c)\), and indeed every scalar-tilt-only quotient
   height, misses part of the \(B\)--\(Q\) synergy.
2. Synergy can be positive with zero angular variation and zero geometric
   projection-height defect.
3. Any long-ray inverse theorem must either control the variation of the
   normalized one-dimensional laws \(\nu_q\), include the likelihood-weighted
   height (0.2), or admit an aligned ray-law-heterogeneity branch.
4. A switching upper bound for (0.2) is a genuinely new quotient-geometric
   statement.  It is not supplied by log-concavity of each ray, by the
   Gaussian channel normalization, or by the finite Bregman identity.

Accordingly, the concrete output is option (C) for the projection-height
mechanism, together with the exact general replacement (2.13), which uses
no quotient expansion.  A proof that
the right side of (2.13) is \(o(L)\) on the positive-mass long-ray family
would close this particular synergy obstruction; no such upper bound is
assumed here.

## 5. Averaging the synergy height does not remove the obstruction

Let \(A:L^2(\mu)\to L^2(\mathbb P_C)\) be the channel conditional-expectation
operator

\[
                         (AF)(c)=\mathbb E[F(X)\mid C=c],
\]

and let \(T=A^*A\) be the posterior-resampling operator on \(L^2(\mu)\).
Thus \(g=AB\) and

\[
                         (TB)(x)=\mathbb E[g(C)\mid X=x]. \tag{5.1}
\]

There is an exact formula for the observation average of (0.2).

**Proposition 5.1 (averaged synergy height).**  For
\(\eta\)-almost every \(q\),

\[
 \boxed{\quad
 \overline h(q):=\mathbb E_C h_C^{\rm syn}(q)
 =L\left\{\frac12-\mathbb E[g(C)\mid Q=q]\right\}
 =L\,\mathbb E[(I-T)B\mid Q=q].
 \quad}                                                   \tag{5.2}
\]

In particular, \(\int\overline h\,d\eta=0\), and

\[
 \left\|{\overline h\over L}\right\|_{L^2(\eta)}^2
 \le\|(I-T)B\|_{L^2(\mu)}^2
 \le\langle B,(I-T)B\rangle_{L^2(\mu)}
 =\mathbb E[g(C)(1-g(C))].                              \tag{5.3}
\]

**Proof.**  The identity \(p_C(c)w_c(q)=p(c\mid q)\) and posterior
calibration give

\[
 \begin{aligned}
 \mathbb E_C[w_C(q)b_q(C)]
 &=\int p(c\mid q)\mathbb P(B=1\mid q,c)dc\\
 &=\mathbb P(B=1\mid Q=q)=\frac12,\\
 \mathbb E_C[w_C(q)g(C)]
 &=\int p(c\mid q)g(c)dc
 =\mathbb E[g(C)\mid Q=q].
 \end{aligned}                                           \tag{5.4}
\]

This proves the first equality in (5.2).  The second follows from (5.1) and
\(\mathbb E[B\mid Q]=1/2\).  Conditional expectation is an \(L^2\)
contraction.  Since \(T=A^*A\) is a positive self-adjoint contraction,
\((I-T)^2\preceq I-T\), proving the two inequalities in (5.3).  Finally,

\[
 \langle B,(I-T)B\rangle
 =\mathbb E B-\mathbb E g(C)^2
 =\mathbb E[g(C)(1-g(C))].
\]

\(\square\)

Formula (5.2) is a first-moment identity and can cancel completely.  In the
diamond model, reflection in the \(T\)-coordinate gives

\[
 g(c_1,-c_2)=1-g(c_1,c_2).                              \tag{5.5}
\]

Conditionally on \(Q=q\), \(C_2=T+G_2\) is symmetric and is independent of
\(C_1=q+G_1\).  Hence

\[
                         \mathbb E[g(C)\mid Q=q]=\frac12
 \quad\hbox{and}\quad
                         \overline h(q)=0               \tag{5.6}
\]

for every \(|q|<R\), despite (3.20).  Jensen averaging of the random
heights therefore erases the entire obstruction in this calibrated
log-concave model.  The common-latent and posterior-resampling identities
see (5.2), not the nonnegative conditional variance \(\zeta\).

There is one likelihood-ratio-free upper bound for the integrated minimal
switching defect, but it has exactly the non-closing scale.  The admissible
correction \(r=h_C^{\rm syn}(Q)\) gives

\[
 \begin{aligned}
 \mathbb E_C\mathcal D_f(h_C^{\rm syn})
 &\le \mathbb E_C\int|h_C^{\rm syn}|d\eta\\
 &=L\,\mathbb E|b_Q(C)-g(C)|\\
 &\le L\sqrt{\mathbb E_C\zeta(C)}.                     \tag{5.7}
 \end{aligned}
\]

If \(f\) is a true extremizer and \(Z=\mathbb E_C\zeta(C)\), (2.12) and
Jensen also give the opposite bound

\[
 {\,\kappa_\beta L\over8}Z^2
 \le\mathbb E_C\mathcal D_f(h_C^{\rm syn})
 \le L\sqrt Z.                                         \tag{5.8}
\]

Thus the normalized channel identities allow switching cost of order \(L\)
when synergy is of constant order.  They do not give the \(o(L)\) upper
bound needed in (2.13).

This order cannot be improved by a scale-free manipulation of the channel
identities.  Let \(\mu_\lambda\) be the image of the diamond law under
\(x\mapsto\lambda x\), let

\[
 f_\lambda(\lambda x)=\lambda f(x),\qquad
 t_\lambda=\lambda^{-2}.                               \tag{5.9}
\]

If \(C_\lambda=t_\lambda X_\lambda+\sqrt{t_\lambda}G\), then
\(\lambda C_\lambda=X+G\).  After the identifications
\(q_\lambda=\lambda q\) and \(c=\lambda c_\lambda\), the quantities
\(w_c,b_q,g,\zeta\), all mutual informations, and all posterior-resampling
probabilities are identical to those of the unit-scale diamond.  The
physical synergy height, however, is

\[
 h_{\lambda,c_\lambda}^{\rm syn}(\lambda q)
 =\lambda h_{1,c}^{\rm syn}(q).                       \tag{5.10}
\]

The minimal correction has the exact scaling law

\[
 \boxed{\quad
 \mathcal D_{f_\lambda,\mu_\lambda}
       (h_{\lambda,c_\lambda}^{\rm syn})
 =\lambda\,
   \mathcal D_{f,\mu}(h_{1,c}^{\rm syn}).
 \quad}                                                 \tag{5.11}
\]

Indeed, \(u_\lambda\) is one-Lipschitz on \(\lambda K\) if and only if
\(u(x)=u_\lambda(\lambda x)/\lambda\) is one-Lipschitz on \(K\), and the
\(L^1\) correction scales by \(\lambda\).

For every \(c\in E\), \(b_q(c)\) is nonconstant.  Since \(w_c>0\) and
\(\int w_c(b_q-g)d\eta=0\), the height \(h_{1,c}^{\rm syn}\) is nonconstant.
The continuous function \(f-h_{1,c}^{\rm syn}(Q)\) cannot be one-Lipschitz:
where the nonconstant smooth height has nonzero derivative, its gradient is
\((-h'(q),1)\), of norm strictly larger than one.  The one-Lipschitz
functions form a closed subset of \(L^1(\mu)\) on the compact diamond, so

\[
 \mathcal D_{f,\mu}(h_{1,c}^{\rm syn})>0
 \qquad(c\in E).                                       \tag{5.12}
\]

Consequently the integrated defect is a positive constant at unit scale and
grows exactly like \(\lambda\) in (5.11), while every dimensionless
common-latent and posterior-resampling statistic stays fixed.  The scaled
model is not isotropic when \(\lambda\ne1\); its purpose is the necessary
affine-scale audit.  It proves that those channel identities alone cannot
yield a likelihood-free \(O(1)\), or more generally \(o(L)\), switching
bound.  Such a bound must use new extremal or quotient geometry.
