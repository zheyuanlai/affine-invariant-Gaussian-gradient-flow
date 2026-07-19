# Gaussian-channel noncoding: exact low-SNR coefficients and the analytic barrier

## 1. Verdict

Let `mu` be isotropic and log-concave, let `S` have mass `1/2`, put

\[
 h=2\mathbf 1_S-1\in\{-1,1\},\qquad
 Y_t=\sqrt tX+G,
 \qquad m_t(Y_t)=\mathbb E[h(X)\mid Y_t],                 \tag{1.1}
\]

and write `g_t=(1+m_t)/2`.  The proposed dimension-free statement

\[
 H(\mathbf1_S\mid Y_{t_0})\geq c                         \tag{1.2}
\]

at a numerical `t_0>0`, uniformly over all such pairs, is not proved here.
Together with `gaussian_profile_localization.md`, a statement of this form
is already a KLS-strength seed.

The low-SNR calculation does give a real improvement over the trace-only
picture.  Define

\[
 a=\mathbb E[hX],\qquad M=\mathbb E[hXX^T].               \tag{1.3}
\]

Then, in nats,

\[
\begin{aligned}
 \chi(t)&:=\mathbb E m_t(Y_t)^2\\
 &=t|a|^2+t^2\left({1\over2}\|M\|_{HS}^2-|a|^2\right)
   +o(t^2),                                                \tag{1.4}\\
 I(h;Y_t)
 &={t\over2}|a|^2+t^2\left({1\over4}\|M\|_{HS}^2
       -{1\over2}|a|^2+{1\over4}|a|^4\right)+o(t^2).      \tag{1.5}
\end{aligned}
\]

The indicator identity and isotropy imply, with exact constants,

\[
             |a|\leq1,\qquad \|M\|_{op}\leq1.             \tag{1.6}
\]

If the dimension-free thin-shell estimate is available for every
log-concave marginal in the squared-radius form

\[
 \operatorname {Var}|Z|^2\leq \tau^2 k                    \tag{1.7}
\]

for every isotropic log-concave `Z in R^k`, then

\[
 \boxed{\ \|M\|_{HS}^2\leq2\tau^2(1+\log n).\ }           \tag{1.8}
\]

Thus the first derivative is dimension free and the second derivative is
only logarithmic, rather than the `O(n)` bound obtained from
`||M||op<=1`.  This is a genuine second-order/M1 gain.

It does **not** produce a seed at `t` of order `1/sqrt(log n)`, or at any
other stated time, because its remainder is not uniform.  The obstruction
is structural: a general log-concave law has only exponential tails.  For a
centered exponential input the Gaussian likelihood ratio is not in
`L^2(gamma)` for any `t>0`, and even its pointwise Taylor series at `t=0`
has radius zero.  Consequently the formal all-chaos Hermite summation which
would propagate (1.4) is invalid.  Log-concave `psi_1` moment estimates do
not repair it.

The stress tests separate this analytic failure from actual coding.

* Every Gaussian partition satisfies the exact no-code estimate
  `chi(t)<=t/(1+t)` by the Ornstein--Uhlenbeck spectrum.  A radial Gaussian
  half-set starts only in the second chaos and has
  `chi(t)=2t^2/pi+o(t^2)`.
* Product cubes also have dimension-free noncoding by tensorization of
  scalar maximal correlation.
* An isotropic atomic codebook `X=+-sqrt(n)e_J` really does reveal its sign
  bit at `t=C log(n)/n`.  It has the right covariance but is not log-concave.
* For the genuine isotropic log-concave product of centered exponentials,
  the balanced max-tail partition retains a numerical amount of conditional
  entropy at **every** sequence `t_n=o(1)`.  This is proved below.  It also
  exhibits the zero-radius Hermite obstruction and shows why a rare-tail
  chaos calculation can falsely resemble a subconstant-time code.
* Replacing an atomic regular simplex by uniform measure on its convex hull
  destroys vertex concentration.  Its max-cell first moment is only
  `Theta(log n/sqrt(n))`; covariance of the atomic vertices alone therefore
  gives the wrong prediction.

The precise remaining problem is a nonperturbative bound on the weighted
quantity `int r_t^2/p_t dgamma`, not another finite-order derivative at zero.

## 2. Exact channel identities

Let `gamma_n` be standard Gaussian measure and set

\[
 L_t(x,y)=\exp\{\sqrt t\langle x,y\rangle-t|x|^2/2\},
\quad p_t(y)=\mathbb E L_t(X,y),
\quad r_t(y)=\mathbb E[h(X)L_t(X,y)].                       \tag{2.1}
\]

The density of `Y_t` relative to `gamma_n` is `p_t`, and Bayes' formula
gives

\[
 m_t(y)={r_t(y)\over p_t(y)},\qquad |r_t(y)|\leq p_t(y).     \tag{2.2}
\]

It follows that

\[
 \boxed{\ \chi(t)=\int {r_t^2\over p_t}\,d\gamma_n.\ }    \tag{2.3}
\]

This denominator is essential for exponential-tailed inputs.  Also,

\[
 \mathbb E[g_t(1-g_t)]={1-\chi(t)\over4}.                   \tag{2.4}
\]

For

\[
 \phi(u)={1\over2}\{(1+u)\log(1+u)+(1-u)\log(1-u)\}
         =\sum_{j\geq1}{u^{2j}\over(2j)(2j-1)},             \tag{2.5}
\]

the mutual information and conditional entropy are

\[
 I(h;Y_t)=\int p_t\phi(r_t/p_t)d\gamma_n,
 \qquad H(h\mid Y_t)=\log2-I(h;Y_t).                        \tag{2.6}
\]

The elementary pointwise bounds

\[
 \phi(u)\leq(\log2)u^2,
 \qquad
 h_2((1+u)/2)\geq(\log2)(1-u^2)                            \tag{2.7}
\]

give

\[
 I(h;Y_t)\leq(\log2)\chi(t),\qquad
 H(h\mid Y_t)\geq(\log2)(1-\chi(t)).                      \tag{2.8}
\]

Thus it is enough to keep `chi` a numerical distance below one.

There is also an exact filtering derivative.  In the equivalent observation
`C_t=tX+B_t`, the posterior magnetization is a martingale and

\[
 \chi'(t)=\mathbb E\left|
       \operatorname {Cov}(h,X\mid C_t)\right|^2.            \tag{2.9}
\]

At zero this is `|a|^2`.  The trace estimate bounds the integrand by a
posterior covariance trace and loses `n`; (1.4) identifies the first tensor
which governs how quickly that loss can appear.

## 3. Hermite calculation through second order

Let `G` be standard Gaussian and write `H_A(G)=G^TAG-Tr A`.  The tensor
Hermite generating identity gives

\[
\begin{aligned}
 p_t(y)&=1+{t\over2}(|y|^2-n)+t^{3/2}p_3(y)+O(t^2),\\
 r_t(y)&=\sqrt t\langle a,y\rangle
       +{t\over2}H_M(y)+t^{3/2}r_3(y)+O(t^2).               \tag{3.1}
\end{aligned}
\]

Here `p_3,r_3` are pure third Gaussian chaoses.  Consequently

\[
\begin{aligned}
 {r_t^2\over p_t}
 ={}&t\langle a,y\rangle^2+t^{3/2}\langle a,y\rangle H_M(y)\\
 &+t^2\left\{{1\over4}H_M(y)^2
       +2\langle a,y\rangle r_3(y)
       -{1\over2}(|y|^2-n)\langle a,y\rangle^2\right\}
       +o(t^2).                                             \tag{3.2}
\end{aligned}
\]

Gaussian-chaos orthogonality kills the degree `1--2` and `1--3` cross
terms.  The two nonzero contractions are

\[
 \mathbb E H_M(G)^2=2\|M\|_{HS}^2,
 \quad {1\over2}\mathbb E[(|G|^2-n)\langle a,G\rangle^2]
      =|a|^2.                                               \tag{3.3}
\]

Equations (3.2)--(3.3) prove (1.4).  Since

\[
 p_t\phi(r_t/p_t)={1\over2}{r_t^2\over p_t}
                  +{1\over12}{r_t^4\over p_t^3}+O(r_t^6/p_t^5)  \tag{3.4}
\]

and `E<a,G>^4=3|a|^4`, (1.5) follows.  In derivative form,

\[
\begin{aligned}
 \chi'(0)&=|a|^2,
 &\chi''(0)&=\|M\|_{HS}^2-2|a|^2,\\
 I'(0)&={1\over2}|a|^2,
 &I''(0)&={1\over2}\|M\|_{HS}^2-|a|^2+{1\over2}|a|^4.    \tag{3.5}
\end{aligned}
\]

For compactly supported laws this follows by direct differentiation.  For
log-concave laws, truncate first and pass to the limit using their finite
sixth moments; this proves the stated one-sided Peano expansions.  It does
not give a dimension-uniform remainder.

To prove (1.6), total covariance for the two conditional laws gives
`I succeq aa^T`, hence `|a|<=1`.  Alternatively,
`|a|=E[h<a/|a|,X>]<=1`.  For every unit `u`,

\[
 |u^TMu|=|\mathbb E[h\langle u,X\rangle^2]|
 \leq\mathbb E\langle u,X\rangle^2=1,                     \tag{3.6}
\]

which proves `||M||op<=1`.

## 4. Thin shell controls the full quadratic tensor up to `sqrt(log n)`

The logarithm in (1.8) comes from a deterministic spectral decomposition,
not from a union bound.

Let `A` be positive semidefinite with eigenvalues
`lambda_1>=...>=lambda_n>=0`, put `lambda_{n+1}=0`, and let `P_k` project
onto the first `k` eigenvectors.  Then

\[
 X^TAX-\operatorname {Tr}A
 =\sum_{k=1}^n(\lambda_k-\lambda_{k+1})
       (|P_kX|^2-k).                                       \tag{4.1}
\]

Every `P_kX` is isotropic and log-concave in its `k`-dimensional range.
Minkowski and (1.7) imply

\[
\begin{aligned}
 \|X^TAX-\operatorname {Tr}A\|_2
 &\leq\tau\sum_k(\lambda_k-\lambda_{k+1})\sqrt k\\
 &=\tau\sum_k\lambda_k(\sqrt k-\sqrt{k-1})\\
 &\leq\tau\sqrt{1+\log n}\,\|A\|_{HS}.                   \tag{4.2}
\end{aligned}
\]

For a general symmetric `A`, apply (4.2) to its positive and negative
parts and use
`||A_+||HS+||A_-||HS<=sqrt(2)||A||HS`.  Duality now gives

\[
\begin{aligned}
 \|M\|_{HS}
 &=\sup_{\substack{A=A^T\\\|A\|_{HS}=1}}
       \mathbb E[h(X)(X^TAX-\operatorname {Tr}A)]\\
 &\leq\sqrt2\,\tau\sqrt{1+\log n},                         \tag{4.3}
\end{aligned}
\]

which is (1.8).

This is the strongest conclusion available from the first two derivatives.
The assertion

\[
 \chi(t)\leq t+C\tau^2t^2\log n                            \tag{4.4}
\]

at a specified noninfinitesimal `t` does **not** follow from (1.4).  It would
require a uniform remainder theorem, which the next section rules out for
the naive Hermite method.

## 5. Why `psi_1` does not give an analytic continuation

Take in one dimension `X=Z-1`, where `Z` has density
`e^{-z}1_{z>=0}`.  The Lebesgue density of `sqrt(t)X+G` has an exponential
right tail:

\[
 f_t(y)\asymp_t e^{-y/\sqrt t}\qquad(y\longrightarrow\infty). \tag{5.1}
\]

Therefore its likelihood ratio `p_t=f_t/varphi` satisfies

\[
 \int p_t^2d\gamma_1=\int {f_t(y)^2\over\varphi(y)}dy
 \gtrsim_t\int^\infty e^{y^2/2-2y/\sqrt t}dy=\infty.       \tag{5.2}
\]

Thus neither Parseval nor an `L^2(gamma)` Hermite summation is available at
any positive time.  Nevertheless (2.3) is finite, because `|r_t|<=p_t`.

The failure is already pointwise at `y=0`:

\[
 p_t(0)=\mathbb E e^{-tX^2/2},\qquad
 {1\over k!}{d^k\over dt^k}p_t(0)\bigg|_{0+}
 ={(-1)^k\mathbb E X^{2k}\over2^kk!}.                       \tag{5.3}
\]

For an exponential variable `E X^{2k}` is asymptotic to `(2k)!`.  Hence
the absolute value of the Taylor coefficient grows like
`(2k)!/(2^kk!)`, whose `k`th root tends to infinity.  The Taylor radius is
zero.  This is fully compatible with the usual log-concave `psi_1` bound
`(E|X|^q)^{1/q}<=Cq`: that bound predicts, rather than prevents, this
growth.

Finite-order Hermite differentiation is legitimate; summing it at a time
depending on `n` is the invalid step.

## 6. Calibration and countertests

### 6.1 Arbitrary Gaussian partitions

For `X~N(0,I_n)`, the normalized observation
`Y_t/sqrt(1+t)` has correlation

\[
 \rho=\sqrt{t/(1+t)}                                      \tag{6.1}
\]

with `X`.  Conditional expectation is the Ornstein--Uhlenbeck operator.
If `h=sum_{k>=1}h_k` is its Gaussian Hermite decomposition, then

\[
 \chi(t)=\sum_{k\geq1}\rho^{2k}\|h_k\|_2^2
 \leq\rho^2={t\over1+t}.                                  \tag{6.2}
\]

This holds for every balanced measurable partition, however oscillatory.
In particular

\[
 H(h\mid Y_t)\geq{\log2\over1+t}.                         \tag{6.3}
\]

For the radial half-set
`h=sign(|X|^2-m_n)`, symmetry gives `a=0` and `M=beta_n I`.  The chi-square
CLT gives

\[
 \beta_n={2\over\sqrt{\pi n}}+o(n^{-1/2}),\qquad
 \|M\|_{HS}^2\longrightarrow{4\over\pi}.                 \tag{6.4}
\]

Hence

\[
 \chi(t)={2\over\pi}t^2+o(t^2),\qquad
 I(h;Y_t)={1\over\pi}t^2+o(t^2).                           \tag{6.5}
\]

Equivalently, the normalized input and output squared radii have limiting
correlation `t/(1+t)`.  The radial cut is a quadratic, not a linear, code.

### 6.2 Product cube

Let the coordinates of `X` be independent uniform variables on
`[-sqrt(3),sqrt(3)]`.  If `T_t` is the scalar conditional-expectation
operator and `rho_cube(t)<1` its norm on the mean-zero scalar subspace, then
the product operator is `T_t^{tensor n}`.  Its norm on the mean-zero product
space is still `rho_cube(t)`, not `sqrt(n)rho_cube(t)`.  Therefore every
balanced cube partition satisfies

\[
 \chi(t)\leq\rho_{cube}(t)^2.                              \tag{6.6}
\]

Boundedness of the scalar input gives `rho_cube(t)->0` as `t->0` (the
kernel converges in Hilbert--Schmidt norm to the independent kernel).
Thus dimension cannot turn a vanishing scalar SNR into a bit on the product
cube.  For the coordinate cut, `a=(sqrt(3)/2)e_1`, `M=0`, so
`chi(t)=3t/4+O(t^2)`.

### 6.3 Covariance-only codebooks really can code

Let `J` be uniform on `{1,...,n}`, let `epsilon` be an independent uniform
sign, and put

\[
 X=\epsilon\sqrt n e_J,qquad h=\epsilon.                   \tag{6.7}
\]

Then `EX=0` and `EXX^T=I`.  Decode by taking the coordinate of largest
absolute value in `Y_t` and then its sign.  With `A=sqrt(tn)`, a union bound
gives

\[
 P\{\widehat h\ne h\}\leq2(n+1)e^{-A^2/8}.                \tag{6.8}
\]

Thus at `t=24 log(n)/n`, the error and `H(h|Y_t)` tend to zero, while
`chi(t)->1`.  The same phenomenon occurs for an atomic isotropic regular
simplex.  These laws are not log-concave; smoothing the separated atoms
produces a Gaussian mixture, still not log-concave.  This is the exact
reason covariance alone cannot prove noncoding.

### 6.4 A rigorous product-exponential max test

Let `Z_i` be iid rate-one exponentials, `X_i=Z_i-1`, and choose `L=L_n` by

\[
 (1-e^{-L})^n={1\over2},\qquad q=e^{-L}=1-2^{-1/n}.          \tag{6.9}
\]

The product law is isotropic and log-concave.  Take

\[
 S_n=\{\max_iZ_i\geq L\}.                                  \tag{6.10}
\]

This is a half-set and `L=log n-log(log2)+o(1)`.  Its perimeter tends to
`(log2)/2`, so it is a genuine calibration, not a KLS counterexample.

**Proposition 6.1.**  For every sequence `t_n->0`, for all sufficiently
large `n`,

\[
 \mathbb E[g_{t_n}(1-g_{t_n})]\geq {1\over2048},
 \quad
 H(\mathbf1_{S_n}\mid Y_{t_n})\geq{\log2\over512}.         \tag{6.11}
\]

**Proof.**  Given the observation, the coordinates remain independent.
The sufficient scalar natural parameter is

\[
 c_i=\sqrt tY_i+t=tZ_i+\sqrt tG_i,                          \tag{6.12}
\]

and the scalar posterior is proportional on the half-line to

\[
 \exp\{(c_i-1)z-tz^2/2\}.                                  \tag{6.13}
\]

Write

\[
 q_i=P\{Z_i\geq L\mid Y_i\},\qquad Q=\sum_iq_i.
\]

Then exactly

\[
 g_t=1-\prod_i(1-q_i),\qquad EQ=nq\leq\log2.               \tag{6.14}
\]

We first show

\[
                  \max_iq_i=o_P(1)                          \tag{6.15}
\]

for every `t=o(1)`.  Put `kappa=tL` and
`u=1+kappa-1/4`.  The moment generating function

\[
 E e^{sc_i}={e^{s^2t/2}\over1-st},\qquad 0<s<1/t,          \tag{6.16}
\]

gives `P(max_i c_i>u)=o(1)`.  Here are explicit Chernoff choices.  In terms
of `r=st`, use `r=1/2` when `kappa<=1/4`, `r=3/4` when
`1/4<=kappa<=1/2`, and

\[
 r=1-{t\over\kappa-1/4}
\]

when `kappa>=1/2`.  After adding `log n=L`, the first two exponents are at
most `-1/(8t)`; in the last regime they are at most

\[
 -{1\over4t}+1+{t\over2(\kappa-1/4)^2}
       +\log{\kappa-1/4\over t},                            \tag{6.17}
\]

which tends to minus infinity.

The posterior tail in (6.13) is increasing in `c`.  If `kappa<=1/4`,
substitution of `c=u<=1` and integration of the decreasing density gives

\[
 q_i\leq C e^{-L/8}.                                       \tag{6.18}
\]

If `kappa>=1/4`, complete the square.  The untruncated Gaussian mode is at
`(c-1)/t`, at least `1/(4t)` to the left of `L`; hence

\[
 q_i\leq C e^{-1/(32t)}.                                   \tag{6.19}
\]

Equations (6.17)--(6.19) prove (6.15).

By Markov's inequality, `P{Q>4log2}<=1/4`.  Combining this with (6.15),
the event

\[
 E_n=\{Q\leq4\log2,\ \max_iq_i\leq1/2\}                  \tag{6.20}
\]

has probability at least `3/4-o(1)`.  On this event,

\[
 1-g_t=\prod_i(1-q_i)\geq e^{-2Q}\geq e^{-8\log2}=2^{-8}. \tag{6.21}
\]

Since `Eg_t=1/2`,

\[
 E[g_t\mathbf1_{E_n}]
 \geq {1\over2}-P(E_n^c)\geq {1\over8}                    \tag{6.22}
\]

for all sufficiently large `n`.  Therefore

\[
 E[g_t(1-g_t)]\geq2^{-8}E[g_t\mathbf1_{E_n}]
 \geq2^{-11}.                                             \tag{6.23}
\]

Finally `h_2(s)>=4(log2)s(1-s)` gives (6.11).  QED.

This example also explains the analytic trap.  Its exact first moment is

\[
 a_i={Lq\over1-q},\qquad |a|^2\sim{(\log2)^2L^2\over n},    \tag{6.24}
\]

and the transverse eigenvalue of `M` is

\[
 {L^2q\over(1-q)^2},                                      \tag{6.25}
\]

so `||M||HS^2=O(L^4/n)`.  The first two chaoses vanish.  Formal high
moments contain factors `qL^k`, however, and an unweighted Hermite sum
suggests a transition near `t=1/L`.  Proposition 6.1 proves that this is not
a code.  The false transition comes from dropping the denominator `p_t` in
(2.3), precisely the operation forbidden by (5.2).

### 6.5 Atomic versus uniform simplex

An atomic regular simplex with `N` vertices is another isotropic codebook:
nearest-vertex decoding succeeds once `tN` dominates `log N`.  It is not a
log-concave probability on its affine span.

For comparison, let `U` be uniform on the simplex
`{u_i>=0, sum u_i=1}` and put

\[
 X=\sqrt{N(N+1)}\left(U-{\mathbf1\over N}\right)           \tag{6.26}
\]

in the `(N-1)`-dimensional zero-sum space.  This law is isotropic and
log-concave.  Let `I=argmax_iU_i`, split the indices into two equal groups,
and let `h` be the group label of `I`.  Since

\[
 E\max_iU_i={H_N\over N},                                  \tag{6.27}
\]

permutation symmetry gives

\[
 \mathbb E[X\mid I=i]
 =\sqrt{N(N+1)}{H_N-1\over N-1}
      \left(e_i-{\mathbf1\over N}\right),                 \tag{6.28}
\]

and hence

\[
 |a|={\sqrt{N+1}(H_N-1)\over N-1}
      =\Theta\left({\log N\over\sqrt N}\right).           \tag{6.29}
\]

Uniform simplex mass is therefore not concentrated at its vertices; the
linear code coefficient tends to zero.  The exponential representation
`U_i=Z_i/sum_jZ_j` makes the connection with Section 6.4 explicit.  No
subconstant-time decoding is proved for this simplex partition, and the
atomic nearest-vertex argument does not transfer to it.

## 7. What is proved and what remains

The rigorous gain is

\[
 \chi'(0)\leq1,
 \qquad
 \chi''(0)\leq2\tau^2(1+\log n),                           \tag{7.1}
\]

with the exact formulas (3.5).  It improves the trace-only infinitesimal
picture and identifies radial/quadratic modes as harmless up to a logarithm.

The rigorous no-go is that no argument based on an `L^2(gamma)` Hermite
radius, or on summing the `psi_1` moments of the likelihood ratio, can turn
(7.1) into a finite-time seed.  Exponential inputs make that radius zero.
The product-exponential max calculation further shows that the apparent
`1/log n` rare-tail transition of the unweighted chaos series is spurious.

An actual continuation needs one of the following genuinely new inputs:

1. a direct inequality for the weighted triangular discrimination
   `int r_t^2/p_t dgamma`, stable under exponential tails;
2. a posterior occupation theorem preventing simultaneous binary certainty
   across many weak directions; or
3. a geometric use of low boundary/perimeter.  Requiring (1.2) for every
   half-set ignores this last information and is stronger than the final
   Cheeger statement.

No universal numerical-time residual uncertainty, and no globally
log-concave subconstant-time code, is established here.
