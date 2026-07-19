# The argmax-sign partition of the isotropic ball does not obstruct projected localization

## 0. Conclusion

Let \(\mu_n\) be the uniform probability measure on

\[
K_n=\sqrt{n+2}\,B_2^n.
\]

This is the isotropic Euclidean ball.  Let \(j(x)\) be the least index at
which \(\max_i|x_i|\) is attained and put

\[
E_n=\{x:x_{j(x)}>0\}.
\tag{0.1}
\]

Ties and the origin have Lebesgue measure zero, so this Borel tie rule has no
effect on any of the identities below.  The involution \(x\mapsto-x\)
interchanges \(E_n\) and its complement almost everywhere.  In particular,
\(\mu_n(E_n)=1/2\).

The proposed stress test does **not** obstruct the bounded-projection
localization argument.  In fact, a deterministic statement stronger than
what the flow needs holds.  Fix \(0<T<1\).  Suppose

\[
d\nu(x)=Z^{-1}\exp\{\langle c,x\rangle-	frac12\langle Bx,x\rangle\}
\mathbf 1_{K_n}(x)\,dx,
\qquad 0\preceq B\preceq TI,
\tag{0.2}
\]

and let \(v\) be any unit eigenvector of \(B\).  Then

\[
\operatorname{Var}_{\nu}\langle X,v\rangle
\le \frac{6}{1-T}.
\tag{0.3}
\]

The bound is uniform in \(n,c,B,v\).  It is proved in Section 3 from one
application of Prekopa--Leindler and one divergence identity.  It does not
use the partition.

For the hard or soft projected localization flow at time \(t\),

\[
B_t=\int_0^t C_s^2\,ds\preceq tI,
\]

and the adaptive survivor \(v_t\), defined as a top eigenvector of
\(K_t=tI-B_t\), is an eigenvector of \(B_t\).  Taking \(T=t=1/2\) in
(0.3) therefore gives the pathwise estimate

\[
\operatorname{Var}_{\mu_{1/2}}\langle X,v_{1/2}\rangle\le 12,
\qquad
\left(2+\operatorname{Var}_{\mu_{1/2}}
\langle X,v_{1/2}\rangle\right)^{-1/2}\ge\frac1{\sqrt{14}}.
\tag{0.4}
\]

For the hard flow the label mass is exactly \(1/2\), on every interval on
which the feedback is well posed.  For the globally well-posed soft flow,

\[
\mathbb P\left\{\frac14\le \mu_s(E_n)\le\frac34
\text{ for every }0\le s\le\frac12\right\}
\ge 1-8\varepsilon,
\tag{0.5}
\]

and the terminal-time probability is at least \(1-2\varepsilon\).  Thus
(0.4) holds surely and it holds simultaneously with balanced survival with
the explicit high probabilities in (0.5).  This proves, rather than merely
suggests, that this model cannot be a high-probability adaptive-survivor
obstruction at the subcritical horizon needed by the mechanism.

There are two superficially alarming but different phenomena.

* Conditioning on the identity of the largest coordinate and then retaining
  both of its signs gives a balanced two-cell mixture with variance
  \(\asymp\log n\) in that coordinate.  This is not the posterior of the
  localization flow at a fixed subcritical time.
* At the critical time \(t=1\), the zero-driving deterministic skeleton with
  a fixed protected diagonal has survivor variance \(\asymp\sqrt n\), while
  its label mass remains exactly balanced.  The exact zero Brownian path has
  probability zero, and the deterministic bound (0.3) deliberately stops at
  every fixed \(T<1\).  Balance alone therefore does not control the critical
  state, but projected localization does not need to run to that state.

Sections 2 and 5 prove these assertions and keep them separate from the
actual pathwise bound (0.4).

---

## 1. The isotropic ball and the measurable label

For the uniform probability measure on the radius-\(R\) Euclidean ball,
rotational symmetry and a one-variable radial integral give

\[
\mathbb E X=0,
\qquad
\mathbb E XX^T=\frac{R^2}{n+2}I.
\]

Consequently \(R=\sqrt{n+2}\) gives the isotropic ball \(K_n\).

Define

\[
j(x)=\min\left\{i:|x_i|=\max_{1\le k\le n}|x_k|\right\}.
\tag{1.1}
\]

The sets in (1.1) are described by finitely many weak inequalities, so
\(j\) is Borel.  Define \(E_n\) by (0.1), declaring the origin not to lie in
\(E_n\).  Write

\[
\sigma(x)=2\mathbf 1_{E_n}(x)-1
\]

outside the null set of ties and coordinate hyperplanes.  Then
\(\sigma(-x)=-\sigma(x)\).  Since \(\mu_n\) is centrally symmetric,

\[
\mu_n(E_n)=\frac12.
\tag{1.2}
\]

Let

\[
L_n=\max_{1\le i\le n}|X_i|.
\]

The initial label-position covariance can be computed exactly.  If
\(j(X)=k\), then \(\sigma(X)X_k=|X_k|=L_n\).  If \(j(X)\ne k\), flipping
only the sign of \(X_k\) preserves the winning coordinate and the label and
cancels the expectation.  Permutation symmetry then gives

\[
\mathbb E[\sigma(X)X_k]
=\mathbb E[L_n\mathbf 1_{\{j(X)=k\}}]
=\frac{\mathbb E L_n}{n}.
\tag{1.3}
\]

Hence, with \({\bf 1}=(1,\ldots,1)^T\),

\[
b_0:=\operatorname{Cov}_{\mu_n}(\mathbf 1_{E_n},X)
=\frac{\mathbb E L_n}{2n}{\bf 1},
\qquad
\frac{b_0}{|b_0|}=\frac{{\bf 1}}{\sqrt n}.
\tag{1.4}
\]

The standard maximum estimates are

\[
c\sqrt{\log(n+1)}\le \mathbb E L_n
\le C\sqrt{\log(n+1)},
\qquad
c\log(n+1)\le \mathbb E L_n^2
\le C\log(n+1).
\tag{1.5}
\]

For completeness, the upper bounds follow from the one-coordinate density

\[
f_n(s)=a_n\left(1-\frac{s^2}{n+2}\right)^{(n-1)/2}
\mathbf 1_{\{|s|\le\sqrt{n+2}\}},
\]

which gives \(\mathbb P\{|X_i|>u\}\le C e^{-cu^2}\), followed by a union
bound and integration of the tail.  For the lower bounds use the exact
spherical representation

\[
X\ \stackrel d=\ \sqrt{n+2}\,
\frac{(G_1,\ldots,G_n)}
{\sqrt{G_1^2+\cdots+G_n^2+H}},
\tag{1.6}
\]

where the \(G_i\) are independent standard Gaussians and \(H\) is an
independent \(\chi^2_2\) variable.  Gaussian and chi-square Chernoff bounds
show, with probability bounded below uniformly for all sufficiently large
\(n\), that the denominator squared in (1.6) is at most \(2(n+2)\) and
\(\max_i|G_i|\ge c\sqrt{\log n}\).  This gives the lower bounds in (1.5);
changing the constants handles the finitely many remaining dimensions.
In particular,

\[
|b_0|\asymp\sqrt{\frac{\log(n+1)}n}.
\tag{1.7}
\]

The small size in (1.7) causes no problem for the soft flow: when the soft
driver initially behaves more like the identity, it adds curvature in the
would-be protected direction and can only move the terminal matrix farther
from the critical rank-one-defect configuration.

---

## 2. Static conditioning is not the adaptive survivor

The even matrix-valued function \(XX^T\) and the odd label \(\sigma\) give

\[
\mathbb E[XX^T\mid E_n]=\mathbb E XX^T=I.
\tag{2.1}
\]

Equation (1.3) gives

\[
m_E:=\mathbb E[X\mid E_n]
=\frac{\mathbb E L_n}{n}{\bf 1}.
\tag{2.2}
\]

Therefore

\[
\operatorname{Cov}(X\mid E_n)
=I-m_Em_E^T\preceq I.
\tag{2.3}
\]

Thus the static half \(E_n\) has no direction with variance of order
\(\log n\).  In particular, conditioning on the entire union in (0.1) must
not be confused with identifying which coordinate wins.

To see exactly where the tempting \(\log n\) scale lives, set

\[
C_{k,+}=\{j(X)=k,\ X_k>0\},
\qquad C_{k,-}=-C_{k,+}.
\]

Each cell has mass \(1/(2n)\).  Conditional on the balanced union
\(C_{k,+}\cup C_{k,-}=\{j(X)=k\}\), the \(k\)-th coordinate has mean zero and

\[
\operatorname{Var}(X_k\mid j(X)=k)
=\mathbb E[L_n^2\mid j(X)=k]
=\mathbb E L_n^2
\asymp\log n.
\tag{2.4}
\]

The second equality follows because \(L_n\) is permutation invariant and
the winning index is uniform.  Formula (2.4) is a rigorous balanced
two-cell \(\log n\) effect.  It is produced by a hard conditioning on a
nonconvex cone union.  It says nothing by itself about the exponentially
quadratically tilted posterior along a projected-localization path.  The
next two sections give a pathwise bound for that posterior.

---

## 3. A subcritical covariance theorem for the ball

We first prove the deterministic analytic lemma used in the conclusion.

### Lemma 3.1 (strong log-concavity of an eigen-direction marginal)

Let \(n\ge2\), \(R>0\), \(0\preceq B\preceq TI\), and let \(v\) be a unit
eigenvector of \(B\), with \(Bv=\beta v\).  Let \(\nu\) have density

\[
d\nu(x)=Z^{-1}e^{\langle c,x\rangle-\langle Bx,x\rangle/2}
\mathbf 1_{\{|x|\le R\}}\,dx.
\tag{3.1}
\]

Put \(d=n-1\).  If

\[
\kappa:=\beta+\frac d{R^2}-T>0,
\tag{3.2}
\]

then the law of \(S=\langle X,v\rangle\) is \(\kappa\)-strongly
log-concave on \((-R,R)\).  Consequently

\[
\operatorname{Var}_{\nu}S\le\frac1\kappa.
\tag{3.3}
\]

#### Proof

Use coordinates \(x=sv+y\), where \(y\in v^\perp\).  Since \(v\) is an
eigenvector, there is no \(s\)-\(y\) cross term.  Write

\[
c=av+h,\qquad B=\beta vv^T+D,
\]

where \(0\preceq D\preceq TI\) on the \(d\)-dimensional space \(v^\perp\).
For \(q>0\), define

\[
G(q)=\int_{|y|^2\le q}
\exp\{\langle h,y\rangle-\tfrac12\langle Dy,y\rangle\}\,dy,
\qquad H(q)=\log G(q).
\tag{3.4}
\]

The set

\[
\{(y,q):|y|^2\le q\}
\]

is convex, and the exponential factor in (3.4) is log-concave in \(y\).
Prekopa--Leindler therefore implies that \(G\) is log-concave as a function
of \(q\).  Equivalently,

\[
H''(q)\le0
\tag{3.5}
\]

in the distributional sense.

We also need a lower bound for \(H'\).  Set \(r=\sqrt q\), let

\[
F(r)=G(r^2),
\]

and let \(\eta_r\) be the probability measure on \(rB_2^d\) proportional
to the integrand in (3.4).  The divergence theorem applied to the vector
field

\[
y\exp\{\langle h,y\rangle-\tfrac12\langle Dy,y\rangle\}
\]

gives the exact identity

\[
r\frac{F'(r)}{F(r)}
=d+\langle h,\mathbb E_{\eta_r}Y\rangle
-\mathbb E_{\eta_r}\langle DY,Y\rangle.
\tag{3.6}
\]

The centered base density
\(e^{-\langle Dy,y\rangle/2}\mathbf 1_{\{|y|\le r\}}\) is even.  Hence its
log-Laplace transform \(\phi_r(h)=\log F(r)\), viewed as a function of
\(h\), is even and convex.  It has its minimum and zero gradient at the
origin, so convexity yields

\[
\langle h,\mathbb E_{\eta_r}Y\rangle
=\langle h,\nabla\phi_r(h)\rangle\ge0.
\tag{3.7}
\]

Moreover \(\langle Dy,y\rangle\le T r^2\) on the integration domain.
Equations (3.6)--(3.7) imply

\[
2qH'(q)=r\frac{F'(r)}{F(r)}\ge d-Tq.
\tag{3.8}
\]

The unnormalized density of \(S\) is

\[
e^{as-\beta s^2/2}G(R^2-s^2).
\]

Put \(L(s)=H(R^2-s^2)\).  From (3.5) and (3.8), at every point of
differentiability and hence distributionally,

\[
\begin{aligned}
L''(s)
&=-2H'(R^2-s^2)+4s^2H''(R^2-s^2)\\
&\le -2H'(R^2-s^2)\\
&\le -\frac d{R^2-s^2}+T\\
&\le -\frac d{R^2}+T.
\end{aligned}
\tag{3.9}
\]

Thus the second distributional derivative of the negative log-density of
\(S\) is at least

\[
\beta+\frac d{R^2}-T=\kappa.
\]

This is precisely \(\kappa\)-strong log-concavity.  The one-dimensional
Brascamp--Lieb inequality gives (3.3).  One may justify the endpoint and
distributional formulation by replacing the ball indicator with a smooth
increasing convex barrier and passing monotonically to the limit; all
constants in (3.8)--(3.9) remain unchanged.  \(\square\)

### Corollary 3.2 (dimension-free subcritical bound)

Let \(R^2=n+2\) and \(0<T<1\).  Under the assumptions of Lemma 3.1,

\[
\operatorname{Var}_{\nu}\langle X,v\rangle
\le\frac6{1-T}.
\tag{3.10}
\]

#### Proof

If

\[
n+2\ge\frac6{1-T},
\]

then

\[
\frac{n-1}{n+2}-T
=1-T-\frac3{n+2}\ge\frac{1-T}{2}.
\]

Lemma 3.1, with \(\beta\ge0\), gives variance at most \(2/(1-T)\).
If instead \(n+2<6/(1-T)\), the random variable
\(\langle X,v\rangle\) is supported on \([-R,R]\).  Popoviciu's inequality
gives variance at most \(R^2=n+2<6/(1-T)\).  This second argument also
covers \(n=1\).  \(\square\)

The threshold \(T=1\) is not an artifact.  Section 5 exhibits a critical
configuration with variance of order \(\sqrt n\).

---

## 4. Application to the hard and soft projected flows

Let \(p_t=d\mu_t/d\mu_n\) be the exponential-quadratic posterior

\[
p_t(x)=Z_t^{-1}\exp\{\langle c_t,x\rangle
-\tfrac12\langle B_tx,x\rangle\}.
\tag{4.1}
\]

For either projected driver,

\[
dc_t=C_t\,dW_t+C_t^2a_t\,dt,
\qquad dB_t=C_t^2\,dt.
\tag{4.2}
\]

For the hard driver, \(C_t=I-u_tu_t^T\) whenever
\(b_t=\operatorname{Cov}_{\mu_t}(\mathbf 1_{E_n},X)\ne0\).  For the soft
driver,

\[
C_t=I-\frac{b_tb_t^T}{|b_t|^2+\varepsilon}.
\tag{4.3}
\]

In both cases \(0\preceq C_t^2\preceq I\), so deterministically

\[
0\preceq B_t\preceq tI.
\tag{4.4}
\]

Write

\[
K_t=tI-B_t.
\]

The survivor \(v_t\) in the rank-one-defect argument is a measurable unit
top eigenvector of \(K_t\).  Since \(B_t=tI-K_t\), the same \(v_t\) is an
eigenvector of \(B_t\).  Corollary 3.2 applies path by path to (4.1).
For every deterministic \(0<t\le T<1\),

\[
\operatorname{Var}_{\mu_t}\langle X,v_t\rangle
\le\frac6{1-T}.
\tag{4.5}
\]

No independence between \(v_t\) and the posterior is used.  This point is
essential: (4.5) directly controls the future-selected direction for which
the generic law-of-total-variance argument fails.

Taking \(t=T=1/2\) gives

\[
\operatorname{Var}_{\mu_{1/2}}\langle X,v_{1/2}\rangle\le12
\tag{4.6}
\]

and hence the pointwise survivor gate

\[
\left(t^{-1}+\operatorname{Var}_{\mu_t}
\langle X,v_t\rangle\right)^{-1/2}
\ge\frac1{\sqrt{14}}
\quad (t=1/2).
\tag{4.7}
\]

For the hard driver, the label mass

\[
M_t=\mu_t(E_n)
\]

satisfies \(dM_t=\langle b_t,C_t\,dW_t\rangle=0\).  Hence

\[
M_t=\frac12
\tag{4.8}
\]

on every interval on which the hard feedback is well posed.  The hard
coefficient is discontinuous at \(b=0\), so (4.8) is not a claim of global
strong well-posedness through a zero.  It remains valid for every stopped
hard solution and for every relaxed zero-signal limit described in the
bounded-projection framework.

For the soft driver, compactness of \(K_n\) makes all posterior moments
globally finite and the smooth parameter SDE is nonexplosive.  The exact
quadratic-variation estimate is

\[
\langle M\rangle_t\le\frac{\varepsilon t}{4}.
\tag{4.9}
\]

Doob's inequality and the terminal Chebyshev bound therefore give

\[
\mathbb P\left\{\sup_{0\le s\le t}|M_s-	frac12|>\tfrac14\right\}
\le16\varepsilon t,
\tag{4.10}
\]

\[
\mathbb P\left\{|M_t-	frac12|>\tfrac14\right\}
\le4\varepsilon t.
\tag{4.11}
\]

At \(t=1/2\), (4.10)--(4.11) become \(8\varepsilon\) and
\(2\varepsilon\), respectively.  Combining them with the sure event (4.6)
proves the advertised high-probability statement (0.5).  In particular,
the adaptive survivor cannot have a \(\log n\), \(\sqrt n\), or any other
unbounded variance at this fixed subcritical horizon, whether or not the
label mass error event occurs.

---

## 5. The critical zero-driving skeleton really is unbounded

This section shows both that the strict inequality \(T<1\) has content and
that balanced mass alone is insufficient at criticality.

Let

\[
v=\frac{{\bf 1}}{\sqrt n},
\qquad c_t=0,
\qquad B_t=tP_{v^\perp},
\quad 0\le t\le1,
\tag{5.1}
\]

and let \(\nu_t\) be the corresponding posterior.  Its density is

\[
d\nu_t(x)\propto
\exp\left\{-\frac t2|x|^2+rac{t}{2n}
\left(\sum_{i=1}^n x_i\right)^2\right\}
\mathbf 1_{K_n}(x)\,dx.
\tag{5.2}
\]

Central symmetry gives \(\nu_t(E_n)=1/2\) and \(a_t=0\).  Permutation
symmetry gives

\[
\operatorname{Cov}_{\nu_t}(\mathbf 1_{E_n},X)\in\mathbb Rv.
\tag{5.3}
\]

The vector in (5.3) is nonzero and points in the positive \(v\)-direction.
Here is a direct proof.  Condition on the coordinate magnitudes
\(a_i=|x_i|\), and let \(j\) be their winning index.  The sign vector
\(\epsilon\in\{-1,1\}^n\) has conditional weight proportional to

\[
\exp\left\{\frac{t}{2n}
\left(\sum_i\epsilon_i a_i\right)^2\right\}.
\]

The Hubbard--Stratonovich identity expresses this sign law as a mixture,
over a symmetric real parameter \(z\), of independent signs with conditional
means \(\tanh(za_i)\).  Consequently

\[
\mathbb E(\epsilon_i\epsilon_j\mid a)
=\mathbb E_z[\tanh(za_i)\tanh(za_j)]\ge0
\quad (i\ne j).
\]

Since \(\sigma(x)=\epsilon_j\),

\[
\mathbb E\left[\sigma(X)\sum_iX_i\mid a\right]
=a_j+\sum_{i\ne j}a_i
\mathbb E(\epsilon_i\epsilon_j\mid a)>0.
\tag{5.4}
\]

Thus the normalized hard feedback direction is exactly \(u_t=v\) along
(5.1).  With the driving increment set equal to zero, (5.1) is a
deterministic skeleton of the hard feedback equations: \(a_t=0\) makes the
drift of \(c_t\) vanish and \(dB_t=P_{v^\perp}dt\).  It is not a
positive-probability Brownian event.

At \(t=1\), put \(S=\langle X,v\rangle\) and \(d=n-1\).  Integrating the
transverse coordinate \(y\in v^\perp\) gives the exact marginal density

\[
f_{n}(s)\propto
\mathbb P\left\{Q_d\le d+3-s^2\right\}
\mathbf 1_{\{|s|\le\sqrt{d+3}\}},
\tag{5.5}
\]

where \(Q_d\) is chi-square with \(d\) degrees of freedom.  Indeed the
transverse integral is

\[
\int_{|y|^2\le n+2-s^2}e^{-|y|^2/2}\,dy
=(2\pi)^{d/2}\mathbb P\{Q_d\le d+3-s^2\}.
\]

The density in (5.5) is even.  There are universal constants \(c,C>0\)
such that

\[
c\sqrt d\le \operatorname{Var}_{\nu_1}S\le C\sqrt d.
\tag{5.6}
\]

We include the short estimate.  For \(s^2>3\), the standard chi-square
lower-tail inequality gives

\[
\mathbb P\{Q_d\le d-(s^2-3)\}
\le \exp\left\{-\frac{(s^2-3)^2}{4d}\right\}
\tag{5.7}
\]

whenever the threshold is nonnegative; it is zero otherwise.  Integration
of (5.7), after the change of variables \(s=d^{1/4}z\), shows that the
normalizing integral in (5.5) is at most \(Cd^{1/4}\) and its unnormalized
second moment is at most \(Cd^{3/4}\).

For the reverse estimates, write

\[
\frac{Q_d-d}{\sqrt{2d}}
=\frac1{\sqrt d}\sum_{i=1}^d\frac{G_i^2-1}{\sqrt2}.
\]

Berry--Esseen applies with a universal constant because
\(\mathbb E|G_i^2-1|^3<\infty\).  Hence, for all sufficiently large \(d\),
the probability in (5.5) is bounded below by a positive universal constant
throughout

\[
\frac12d^{1/4}\le |s|\le d^{1/4}.
\]

This interval gives a lower bound \(cd^{1/4}\) for the normalizer and
\(cd^{3/4}\) for the unnormalized second moment.  The same upper tail bound
gives upper bounds \(Cd^{1/4}\) and \(Cd^{3/4}\), respectively.  Dividing
the lower second-moment bound by the upper normalizer gives the lower bound
in (5.6), while dividing the upper second-moment bound by the lower
normalizer gives the upper bound.  Enlarging the constants handles bounded
\(d\).  This proves (5.6).

Thus a perfectly balanced, self-consistent deterministic skeleton has an
unbounded adaptive survivor at the critical time.  It does not contradict
Section 4: the projected-localization proof is free to stop at \(t=1/2\),
where (4.6) holds on every stochastic path.  Nor does (5.6) show that the
actual hard process has an unbounded survivor with nonnegligible
probability at \(t=1\); proving such a statement would require a stochastic
small-ball estimate for the driving signal, which (5.1) does not provide.

---

## 6. Small dimensions and asymptotic summary

* **\(n=1\).**  The measure is uniform on
  \([-\sqrt3,\sqrt3]\), \(E_1=(0,\sqrt3]\), and the hard projected driver is
  zero because the protected subspace is the whole line.  The posterior is
  static and the survivor variance is exactly \(1\).  Every soft posterior
  is supported in the same interval and has variance at most \(3\).
* **\(2\le n<6/(1-T)-2\).**  The support estimate
  \(\operatorname{Var}\langle X,v\rangle\le n+2<6/(1-T)\) is already
  dimension-free because only finitely many dimensions occur for fixed
  \(T<1\).
* **Large \(n\).**  Lemma 3.1 gives the sharper estimate
  \(\operatorname{Var}\langle X,v\rangle\le
  (1-T-3/(n+2))^{-1}\), so the bound tends to \((1-T)^{-1}\).
* **Criticality.**  At \(T=1\), the transverse Gaussian radius has
  chi-square fluctuations of order \(\sqrt n\).  The unpenalized chord
  consequently has width \(n^{1/4}\), producing the \(\sqrt n\) variance
  in (5.6).  This is exactly why the denominator in the subcritical bound
  is \(1-T\).

The final distinction is therefore exact.  Static conditioning on the
whole half \(E_n\) has covariance at most the identity.  A refined hard
cell-pair has a \(\log n\) variance.  A critical zero-noise skeleton has a
\(\sqrt n\) variance.  The actual soft or hard projected posterior at the
fixed subcritical time \(1/2\), with its future-selected survivor direction,
has variance at most \(12\) on every path, and the label is balanced with
the explicit probabilities (0.5).  The isotropic-ball argmax-sign model is
therefore rigorously excluded as an obstruction to the bounded-projection
localization mechanism.
