# Continuous-bootstrap and cut-locus audit

This note asks whether the known stochastic-localization bootstraps can be
read as a dimension-free fixed-point argument.  The answer is negative for
the known inequalities.  There are two distinct obstructions:

1. the ordinary covariance process has a genuine `log n` extreme-coordinate
   time scale even when the initial Poincare constant is universal; and
2. the Chen exponent bootstrap has a dimension-dependent positive equilibrium
   when written as a continuous recurrence.

The last section audits a proposed full cut-locus second-variation term.  The
long-ray estimates control smooth curvature on the long core but give no upper
bound on the singular/medial energy, which is measured with a different
weight.  An explicit ray-data model shows that the latter can pay the entire
normal-variance inequality without changing the bad first-moment scale.

Throughout, `C,c` are numerical constants.  Put

\[
 K_n=\sup\{C_P(\mu):\mu\text{ isotropic and log-concave in dimension at
 most }n\}.
\tag{1}
\]

The use of dimensions at most `n` only makes `K_n` monotone and does not
change any conclusion.

## 1. The exact improved-Lichnerowicz fixed point

Let `p` be isotropic and log-concave, and run ordinary Eldan localization

\[
 p_{t,\theta}(x)=Z(t,\theta)^{-1}
 e^{\langle\theta,x\rangle-t|x|^2/2}p(x).
\tag{2}
\]

The random tilt has the filtering representation

\[
 \theta_t\stackrel{d}=tX+W_t,
\tag{3}
\]

where `X~p`, `W_t~N(0,tI)`, and they are independent.  Write `A_t` for
the covariance of `p_(t,theta_t)`.

We use two precisely stated standard inputs.

* **Covariance-survival lemma.**  If

  \[
  \kappa_n^2=\sup_{\mu,u}
  \left\|\mathbb E_\mu[\langle X,u\rangle X\otimes X]\right\|_{HS}^2,
  \tag{4}
  \]

  where the supremum is over isotropic log-concave laws of dimension at
  most `n` and unit `u`, then

  \[
  \mathbb E\|A_t\|_{op}\le C_A
  \quad\text{for }0\le t\le {c_A\over\kappa_n^2\log(en)}.
  \tag{5}
  \]

  This is the covariance estimate used in the stochastic-localization proof
  of the logarithmic KLS bound.  Its hypotheses include no KLS conclusion.

* **Improved log-concave Lichnerowicz inequality.**  If `nu` is
  `t`-uniformly log-concave and has covariance `A`, then

  \[
  C_P(\nu)\le\sqrt{\|A\|_{op}/t}.
  \tag{6}
  \]

For completeness, the self-consistent estimate for (4) follows directly
from Poincare.  Fix `u`, put

\[
 M=\mathbb E[\langle X,u\rangle X\otimes X],
 \qquad F(x)=\langle Mx,x\rangle.
\]

Since `E<X,u>=0`, isotropy and Cauchy--Schwarz give

\[
 \|M\|_{HS}^2
 =\mathbb E[\langle X,u\rangle(F(X)-\mathbb EF(X))]
 \le\sqrt{\operatorname{Var}F(X)}.
\]

Moreover, `grad F=2Mx`, so

\[
 \operatorname{Var}F(X)
 \le C_P(\mu)\mathbb E|2MX|^2
 =4C_P(\mu)\|M\|_{HS}^2.
\]

Cancellation of the nonzero factor (the zero case is trivial) proves

\[
 \boxed{\kappa_n^2\le4K_n.}
\tag{7}
\]

We also spell out the endpoint transfer.  Let `lambda_0=1/C_P(p)`.  The
variance decomposition along (2), the Poincare inequality for the tilt
`theta_t=tX+W_t`, and E. Milman's existence of a 1-Lipschitz function with
variance at least `c_M C_P(p)` give

\[
 C_P(p)\le C(2+tC_P(p))\,\mathbb E C_P(p_t).
\tag{8}
\]

Here is the short derivation.  If `f` is the Milman function and
`M_t=int f p_t`, then

\[
 \operatorname{Var}_p f
 =\mathbb E\operatorname{Var}_{p_t}f+\operatorname{Var}M_t.
\]

The tilt law has Poincare constant at most
`t^2 C_P(p)+t`.  Differentiation in `theta` and covariance
Cauchy--Schwarz yield

\[
 |\nabla_\theta M_f(t,\theta)|^2
 \le\|A(t,\theta)\|_{op}\operatorname{Var}_{p_(t,theta)}f
 \le t^{-1}\operatorname{Var}_{p_(t,theta)}f.
\]

Thus `Var M_t <= (1+t C_P(p)) E Var_(p_t)f`, proving (8).  Combining
(6), (8), and Jensen shows that whenever `t C_P(p)<=a`,

\[
 C_P(p)\le {C_a\over\sqrt t}
       \sqrt{\mathbb E\|A_t\|_{op}}.
\tag{9}
\]

Choose

\[
 t={c\over K_n\log(en)}.
\tag{10}
\]

Equations (5), (7), and `C_P(p)<=K_n` verify all conditions in (9), and
give

\[
 C_P(p)\le C\sqrt{K_n\log(en)}.
\]

Taking the supremum and solving the scalar inequality gives the exact known
fixed point

\[
 \boxed{K_n\le C\sqrt{K_n\log(en)}
        \quad\Longrightarrow\quad K_n\le C^2\log(en).}
\tag{11}
\]

There is already no iteration depth in (11).  Iterating the right-hand side
does not improve its fixed point.  Indeed, if

\[
 K_{j+1}=C\sqrt{K_j L},\qquad L=\log(en),
\]

then

\[
 K_j=C^{2(1-2^{-j})}L^{1-2^{-j}}K_0^{2^{-j}},
\tag{12}
\]

whose limit is `C^2 L`.

If (10) could be improved by replacing `log(en)` with
`1+log K_n`, then (11) would become

\[
 K_n\le C\sqrt{K_n(1+\log K_n)},
\tag{13}
\]

and hence `K_n<=C^2(1+log K_n)`, which forces a numerical upper bound.
Section 3 proves that this proposed covariance improvement is false even for
a product law with universal Poincare constant.

## 2. What dimension-free thin shell supplies, and where it stops

For an isotropic log-concave vector in dimension `k`, put

\[
 S_k=\sup {1\over k}\operatorname{Var}|X|^2.
\tag{14}
\]

The dimension-free thin-shell theorem says `sup_k S_k=S<infinity`.
It implies only the harmonic estimate

\[
 \boxed{\kappa_n^2\le C\sum_{k=1}^n{S_k\over k}
                  \le CS\log(en).}
\tag{15}
\]

Here is the algebra.  For `M` as above, diagonalize it.  If `P` projects
onto any `r` eigenvectors, then the marginal `PX` is isotropic and
log-concave on its range, and

\[
 |\operatorname{Tr}(PM)|
 =|\mathbb E\langle X,u\rangle(|PX|^2-r)|
 \le\sqrt{rS_r}.
\tag{16}
\]

Apply (16) separately to the positive and negative eigenvalues.  The `r`th
eigenvalue in either sign class has square at most `S_r/r`.  Summing both
sign classes proves (15), up to a numerical factor.

Combining (15), (5), and (9) gives

\[
 t_0\asymp {1\over S(\log(en))^2},
 \qquad K_n\le C\sqrt S\log(en).
\tag{17}
\]

Using the classical rather than improved Lichnerowicz transfer gives the
older `K_n<=CS(log(en))^2`.  Combining (7) and (15) merely replaces
`kappa_n^2` by

\[
 \min\{4K_n,CS\log(en)\};
\]

the resulting inequality

\[
 K_n\le C\sqrt{\min\{K_n,\log(en)\}\log(en)}
\tag{18}
\]

still has order `log n` as its smallest uniform fixed point.

Projection thin-shell information cannot remove the harmonic sum in (15).
Let

\[
 B=\operatorname{diag}(1,2^{-1/2},\ldots,n^{-1/2}),
 \qquad \mathcal Q(H)=(\operatorname{Tr}BH)^2.
\tag{19}
\]

For every rank-`r` orthogonal projection `P`, Ky Fan's inequality gives

\[
 \mathcal Q(P)
 \le\left(\sum_{i=1}^r i^{-1/2}\right)^2\le4r.
\tag{20}
\]

On the other hand,

\[
 \sup_{\|H\|_{HS}=1}\mathcal Q(H)
 =\|B\|_{HS}^2=\sum_{i=1}^n{1\over i}\ge\log(n+1).
\tag{21}
\]

Thus all projection inequalities used in (16) are compatible with a
`log n` general quadratic-form norm.  Any improvement of (15) must use
additional compatibility of the third moment tensor, not the thin-shell
theorem by itself.

## 3. An actual log-concave process showing that `log n` is sharp

Let `Y=E-1`, where `E` is a unit exponential variable, and set

\[
 \mu_n=\operatorname{Law}(Y)^{\otimes n}.
\tag{22}
\]

This law is isotropic and log-concave.  The one-dimensional exponential
Poincare constant is at most `4`, and tensorization gives

\[
 C_P(\mu_n)\le4
\tag{23}
\]

for every `n`.

At localization time `t`, the one-dimensional posterior with tilt `theta`
has density

\[
 p_{t,\theta}(y)\propto
 \exp\left((\theta-1)y-{t\over2}y^2\right)
 \mathbf1_{\{y\ge-1\}}.
\tag{24}
\]

Put `m=(theta-1)/t`.  Under (24),

\[
 Z=\sqrt t\,(Y-m)
\]

is standard Gaussian conditioned on

\[
 Z\ge a=-\sqrt t-{\theta-1\over\sqrt t}.
\tag{25}
\]

If `theta>=1+sqrt(t)`, then `a<=-1`.  Let `G` be standard Gaussian and
write

\[
 p_I=P(-1\le G\le0),\qquad p_J=P(1\le G\le2).
\]

For every `a<=-1`, the conditional probabilities of these two intervals
are at least `p_I,p_J`.  The independent-copy formula for variance therefore
gives

\[
 \operatorname{Var}(G\mid G\ge a)\ge p_Ip_J=:c_0>0.
\]

Consequently

\[
 \operatorname{Var}_{p_(t,theta)}Y\ge {c_0\over t}
 \quad\text{when }\theta\ge1+\sqrt t.
\tag{26}
\]

Now use the filtering representation (3).  Coordinates are independent and

\[
 \theta_{t,i}=tY_i+W_{t,i}.
\]

The event

\[
 \{Y_i\ge0,\ W_{t,i}\ge1+\sqrt t\}
\tag{27}
\]

has probability

\[
 p_t=e^{-1}\overline\Phi(t^{-1/2}+1)
\tag{28}
\]

and implies the hypothesis of (26).  Take `t=1/log n`.  The elementary
Mills lower bound

\[
 \overline\Phi(u)\ge {u\over1+u^2}{e^{-u^2/2}\over\sqrt{2\pi}}
\]

shows

\[
 np_t\ge {c\over\sqrt{\log n}}
 \exp\left({\log n\over2}-\sqrt{\log n}\right)\longrightarrow\infty.
\tag{29}
\]

Hence, for all sufficiently large `n`, with probability at least `1/2` at
least one coordinate satisfies (27).  Since the posterior remains a product,
its covariance is diagonal.  Equations (26)--(29) imply

\[
 \boxed{
 \mathbb E\|A_{1/\log n}\|_{op}\ge c\log n,
 \qquad
 \mathbb E\sqrt{\|A_{1/\log n}\|_{op}}\ge c\sqrt{\log n}.}
\tag{30}
\]

This is an obstruction inside the class of isotropic log-concave measures,
not merely a matrix martingale analogy.  In particular, no estimate of the
form

\[
 \mathbb E\|A_t\|_{op}\le C
 \quad(0\le t\le c/[K_n(1+\log K_n)^a])
\tag{31}
\]

can hold: (23) makes the proposed horizon numerical, while (30) fails at a
time tending to zero.  The `log n` in ordinary operator-norm covariance
survival is therefore sharp even on a product for which KLS is elementary.
A successful proof must avoid charging an irrelevant extreme posterior
coordinate to every test function or cut.

At time zero this example has the exact symmetric third tensor

\[
 \mathbb E[Y_iY_jY_k]=2\mathbf1_{\{i=j=k\}}.
\]

Thus the martingale part of the covariance SDE begins as `n` independent
diagonal noises `2E_(ii)dW_i`.  This gives the elementary process-level
reason for (30): the maximum of `n` scalar noises has a `sqrt(log n)`
excursion, and its variance scale is `log n`.

## 4. Chen's bootstrap as a continuous recurrence

Write `R_d=1/psi_d` for the worst normalized inverse Cheeger constant in
dimension at most `d`.  Chen's induction lemma has the following exact
form.  If

\[
 R_k\le\alpha k^\beta\quad(1\le k\le d),
 \qquad0<\beta\le1/2,
\]

and `q=ceil(1/beta)+1`, then

\[
 R_d\le c\,q^{1/2}\alpha(\log d)^{1/2}
          d^{\beta-\beta/(8q)}.
\tag{32}
\]

Since

\[
 {1\over\beta}<q\le {2\over\beta},
 \qquad {\beta^2\over16}\le {\beta\over8q}< {\beta^2\over8},
\tag{33}
\]

the exponent update is

\[
 \beta\longmapsto\beta-\Theta(\beta^2),
\tag{34}
\]

while the prefactor is multiplied by at least the displayed
`sqrt(log d)` and by order `beta^(-1/2)`.  This is why a depth growing with
dimension appears in that argument.

Even grant the bootstrap its most favorable self-consistent
reinterpretation.  Put `L=log d`, `r=log R_d`, take
`beta=r/L`, and suppress all previous prefactor losses.  Equations
(32)--(33) still give only

\[
 r_{new}\le r-{r^2\over16L}+\log L+C
\tag{35}
\]

for `1<=r<=L/2`.  The decrement beats the cost only when

\[
 r\gtrsim\sqrt{L\log L}.
\]

Thus the idealized continuous recurrence has a positive equilibrium

\[
 r_*\asymp\sqrt{\log d\,\log\log d},
\tag{36}
\]

rather than a numerical fixed point.  The rigorous discrete recursion
`beta_(j+1)=beta_j-beta_j^2/16` and
`alpha_(j+1)=2c alpha_j beta_j^(-1/2)` is the same phenomenon without the
favorable suppressions.  No rearrangement of (32) produces an inequality
`R<=C(1+log R)^a`; one would have to replace either the `log d` charge or
the `r^2/log d` gain by a new estimate.

## 5. Audit of a full cut-locus/medial-energy second variation

Let `A=D_1(mu)` and let a true T3 extremizer be disintegrated into calibrated
balanced rays.  Existing extremal structure gives numerical
`alpha_0,c_0>0` and a quotient set `Omega` with

\[
 \eta(\Omega)\ge\alpha_0,
 \qquad \sigma_y\ge c_0A\quad(y\in\Omega).
\tag{37}
\]

On a smooth two-sided ray chart, the focal estimate gives

\[
 \|II_y\|_{HS}^2\le {C\over A^2}
 \quad(y\in\Omega).
\tag{38}
\]

Moreover, the small-event barycenter lemma shows that every coherent unit
cap of long-ray directions has quotient mass at most `C exp(-cA)`.  Hence a
connected completion of the Gauss image must carry exponentially many
unit-scale direction packets, and therefore potentially exponentially large
absolute singular turning.  This strengthens, rather than upper-bounds, the
medial term.

Suppose hypothetically that a rigorous full second-variation formula gave

\[
 \operatorname{Var}_{\widehat\sigma}(N)
 \le C\big(\mathcal E_{smooth}+\mathcal E_{medial}\big),
\tag{39}
\]

where `widehat(sigma)` is normalized boundary measure and the second term
contains the cut-locus graph.  Equation (39) has the wrong direction for a
bootstrap.  The dispersed directions give a lower bound on its right-hand
side.  Equations (37)--(38) upper-bound only the smooth contribution on the
long core; they give no upper bound at all on the complementary or medial
contribution.

The mismatch is not cosmetic.  Here is explicit abstract ray data, satisfying
all presently used one-dimensional, covariance, thin-shell, perimeter-scale,
and long-core focal constraints.  Let `A>=10`, put

\[
 m=\lceil A^2\rceil,\qquad n\ge A^4.
\]

Give total quotient mass `1-A^(-1)` to `m` equally weighted directions
`e_1,...,e_m`.  On each such ray take the uniform law on
`[-sqrt(3)A,sqrt(3)A]`.  Then

\[
 \sigma_y=A,\quad q_y(0)={1\over2\sqrt3 A},\quad II_y=0,
\]

and

\[
 \int\sigma_y^2N_yN_y^T\,d\eta\preceq I,
 \qquad \int\sigma_y^4d\eta\le n.
\tag{40}
\]

Give the remaining quotient mass `A^(-1)` to unit-scale transition rays.
Their boundary density is numerical, so their contribution to the unnormalized
boundary measure is order `A^(-1)`, the same order as the entire long core.
Consequently they carry a numerical fraction of normalized boundary measure,
although they contribute only `O(A^(-1))` to the first-moment objective.
Join the `m` long-ray normals through an arbitrary medial graph on this
transition part.  Its singular energy can be made as large as desired, while
(40), the `D_1` scale `asymp A`, and total perimeter `asymp A^(-1)` do not
change.  In particular it can pay (39) completely.

This is a measure-theoretic version of the long disconnected-cylinder charts
glued through a short fan.  It is not asserted to be globally log-concave;
its purpose is exact: (39), (37), (38), covariance, and thin shell do not
algebraically imply any recurrence for `A`.

A closing estimate would need genuinely new global content of one of the
following forms:

\[
 \mathcal E_{medial}\le {C\over A^p}\quad(p>0),
\tag{41}
\]

or a variational inequality converting medial charge into a competing
1-Lipschitz function whose first moment is strictly larger than `A`.  If
(41) held also for the complementary smooth energy, then a numerical lower
bound on `Var(N)` and (39) would immediately force `A=O(1)`.  No known
second-variation, thin-shell, or localization identity supplies (41).  The
Gaussian fan demonstrates why the singular charge cannot simply be omitted;
the weight model above demonstrates why long-ray flatness does not control it.

## 6. Conclusion

The known localization mechanisms have already been put in their strongest
self-consistent fixed-point form:

\[
 K\mapsto C\sqrt{K\log n}.
\]

Its `O(log n)` fixed point is sharp for operator-norm covariance survival,
as witnessed by the product exponential process.  Thin shell leaves both a
harmonic quadratic-form obstruction and the same extreme-coordinate entropy.
Chen's continuous exponent recursion has equilibrium
`exp(Theta(sqrt(log n log log n)))`.  Finally, a full cut-locus term helps
only if accompanied by a new *upper* bound or a competitor construction for
the medial energy.  None of these missing estimates is a consequence of the
known bootstrap inequalities.
