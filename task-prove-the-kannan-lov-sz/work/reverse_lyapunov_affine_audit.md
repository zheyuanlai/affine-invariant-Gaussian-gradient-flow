# Reverse convolution Lyapunovs and the Gaussian-channel affine gate

## 0. Conclusions

Let

\[
 \mathsf T\mu=\mathcal L((X_1+X_2)/\sqrt2),\qquad
 a=\lambda_1(\mu),\quad b=\lambda_1(\mathsf T\mu).
\]

The forward estimate in the companion note gives
\(b(1-b)\leq4(b-a)\).  A reverse Lyapunov would have to pay for the
total amplification

\[
 \sum_j\log(a_{j+1}/a_j)=-\log a_0.
\]

This audit gives four exact conclusions.

1. Raw relative Fisher information has no universal budget and is infinite
   on basic hard-support models.  A bounded inverse-Fisher deficit is
   destroyed by a Gaussian regularization of the one-sided exponential:
   its decrement tends to zero while the gap amplification tends to
   \(\log(16/9)\).
2. Fixed-noise Gaussian-channel Fisher information, equivalently the MMSE
   deficit, does have a universal budget.  Smooth strongly log-concave
   high-Hermite perturbations nevertheless make the ratio of spectral
   amplification to its decrement grow as \(2^m/(m-2)\).
3. Entropy deficit and Hadamard sum/difference mutual information survive
   both one-dimensional tests locally.  Their scalar budgets are additive
   under products.  Division by dimension fails under Gaussian padding,
   whereas bounded transforms of the raw additive quantities saturate.
4. In the proposed Gaussian-channel proof of mean-gradient coercivity, all
   non-affine posterior modes are controlled by the Hessian.  The exact
   remaining term is

   \[
       \mathbb E[(I-A_Y)a_Y],\qquad
       A_y=\operatorname {Cov}(X\mid Y=y),\quad
       a_y=\mathbb E[\nabla g(X)\mid Y=y].
   \]

   Posterior strong log-concavity does not control this affine term.  The
   desired estimate for it is the load-bearing new input.

All Hermite asymptotics below mean that \(m\) is first fixed and
\(\varepsilon\downarrow0\); only afterwards may \(m\) tend to infinity.

## 1. The standardized Gamma calibration

Let \(A_k\sim\operatorname {Gamma}(k,1)\) and

\[
 X_k={A_k-k\over\sqrt k}.
\]

Then \(X_k\) is centered with variance one, and normalized
self-convolution sends \(X_k\) exactly to \(X_{2k}\).  Its Euclidean
Poincare gap is

\[
 g_k={k^2\over(k+1)^2}.
\]

Consequently one step and the full dyadic amplification are

\[
 \ell_k:=\log{g_{2k}\over g_k}
 =2\log {2(k+1)\over2k+1},\qquad
 \sum_{j\geq0}\ell_{2^jk}=2\log(1+1/k).
 \tag{1.1}
\]

### 1.1 Entropy

Writing \(h_G=\frac12\log(2\pi e)\), the exact entropy and deficit are

\[
 h_k=k+\log\Gamma(k)+(1-k)\psi(k)-{1\over2}\log k,
 \qquad \delta_k=h_G-h_k.
 \tag{1.2}
\]

Stirling expansion gives

\[
 \delta_k={1\over3k}+{1\over12k^2}+{1\over90k^3}
 +O(k^{-4}),\qquad
 \delta_k-\delta_{2k}={1\over6k}+O(k^{-2}).
 \tag{1.3}
\]

Thus \(\ell_k/(\delta_k-\delta_{2k})\to6\).  For \(k=1\),
\(h_1=1\) and \(h_2-h_1=\gamma_{\rm E}-\frac12\log2\).

### 1.2 Unsmooth Fisher information

Inside its support, the score of \(X_k\) is

\[
 \rho_k(X_k)=\sqrt k\left({k-1\over A_k}-1\right).
\]

For \(k>2\),

\[
 \mathbb EA_k^{-1}={1\over k-1},\qquad
 \mathbb EA_k^{-2}={1\over(k-1)(k-2)},
\]

and hence

\[
 J(X_k)={k\over k-2},\qquad
 I(X_k\Vert\gamma)=J(X_k)-1={2\over k-2},
 \tag{1.4}
\]

\[
 I(X_k\Vert\gamma)-I(X_{2k}\Vert\gamma)
 ={k\over(k-2)(k-1)}.
 \tag{1.5}
\]

The Sobolev Fisher information is infinite for \(k\leq2\).  At \(k=1\)
the density has a nonzero boundary trace, and at \(k=2\) it vanishes
linearly, producing \(\int_0dt/t\).  Integrating only the pointwise
interior score at \(k=1\) misses the distributional boundary derivative.

The bounded harmonic deficit \(1-1/J\) equals \(2/k\) for \(k>2\), and
its one-step decrement is \(1/k\).  Section 2 shows why this attractive
calibration cannot survive hard support.

### 1.3 Exact Gaussian-channel MMSE formula

Let \(G\sim N(0,1)\) be independent and \(Y_k=X_k+G\).  In terms of the
parabolic-cylinder function \(D_\nu\), the output density and score are

\[
 q_k(y)={k^{k/2}\over\sqrt{2\pi}}
 \exp\left(-{y^2\over4}-\sqrt k\,y-{k\over2}\right)D_{-k}(-y),
 \tag{1.6}
\]

\[
 \rho_{q_k}(y)={D_{1-k}(-y)\over D_{-k}(-y)}-\sqrt k.
 \tag{1.7}
\]

These follow by completing the square and using

\[
 \int_0^\infty u^{k-1}e^{-u^2/2+yu}du
 =\Gamma(k)e^{y^2/4}D_{-k}(-y).
\]

Define the bounded Gaussian-channel Fisher/MMSE deficit

\[
 \mathcal F_k
 :=I\left({X_k+G\over\sqrt2}\middle\Vert\gamma\right)
 =2\int \rho_{q_k}^2q_k-1
 =1-2\,\operatorname {mmse}(X_k\mid X_k+G).
 \tag{1.8}
\]

The last identity uses
\(\operatorname {mmse}(X\mid X+G)=1-J(X+G)\).
An explicit posterior formula is

\[
 \mathbb E[A_k^r\mid Y_k=y]
 =k^{r/2}(k)_r{D_{-k-r}(-y)\over D_{-k}(-y)}.
 \tag{1.9}
\]

A Hermite calculation gives

\[
 \mathcal F_k={1\over4k}-{1\over8k^2}+O(k^{-3}),\qquad
 \mathcal F_k-\mathcal F_{2k}={1\over8k}+O(k^{-2}).
 \tag{1.10}
\]

Thus \(\ell_k/(\mathcal F_k-\mathcal F_{2k})\to8\).  Direct quadrature of
(1.6)-(1.8) gives

\[
 \mathcal F_1=0.1764315\ldots,\quad
 \mathcal F_2=0.1022906\ldots,\quad
 \mathcal F_4=0.0559862\ldots.
\]

### 1.4 Sum/difference dependence

Let

\[
 S_k={X_k+X_k'\over\sqrt2}=X_{2k},\qquad
 D_k={X_k-X_k'\over\sqrt2}.
\]

The correct mutual-information sign is

\[
 I(S_k;D_k)=h(S_k)+h(D_k)-2h(X_k)\geq0.
 \tag{1.11}
\]

For \(U_k=A_k-A_k'\),

\[
 p_{U_k}(u)=
 {|u|^{k-1/2}K_{k-1/2}(|u|)
  \over2^{k-1/2}\sqrt\pi\,\Gamma(k)},\qquad
 D_k={U_k\over\sqrt{2k}}.
 \tag{1.12}
\]

If \(\delta_{D,k}=h_G-h(D_k)\), then

\[
 I(S_k;D_k)=2\delta_k-\delta_{2k}-\delta_{D,k}.
 \tag{1.13}
\]

The difference has fourth cumulant \(3/k\), so the standard Edgeworth
calculation gives

\[
 \delta_{D,k}={3\over16k^2}+O(k^{-3}),\qquad
 I(S_k;D_k)={1\over2k}-{1\over24k^2}+O(k^{-3}).
 \tag{1.14}
\]

For \(k=1\), \(D_1\) is variance-one Laplace and

\[
 h(X_1)=1,\quad
 h(S_1)=1+\gamma_{\rm E}-{1\over2}\log2,\quad
 h(D_1)=1+{1\over2}\log2.
\]

Therefore

\[
 \boxed{I(S_1;D_1)=\gamma_{\rm E}.}
 \tag{1.15}
\]

The joint pair still has Poincare constant \(4\), although its two visible
marginals have constants \(9/4\) and \(2\).  Dependence carries the
original slow mode.

The axial marginal of the isotropic simplex converges to \(X_1\), and its
\(r\)-fold normalized-sum sector converges to \(X_r\).  Every formula above
is therefore also a mandatory limiting test for the simplex axial sector.

## 2. Hard support destroys bounded inverse-Fisher Lyapunovs

Let \(X_1=E-1\), where \(E\sim\operatorname {Exp}(1)\), and set

\[
 Z_\varepsilon={X_1+\sqrt\varepsilon G\over\sqrt{1+\varepsilon}},
 \qquad \mu_\varepsilon=\mathcal L(Z_\varepsilon).
 \tag{2.1}
\]

This is smooth, full-support, isotropic, and log-concave.  Moreover,

\[
 \mathsf T\mu_\varepsilon
 =\mathcal L\left({X_2+\sqrt\varepsilon G'
                  \over\sqrt{1+\varepsilon}}\right).
 \tag{2.2}
\]

The gaps satisfy

\[
 \lambda_1(\mu_\varepsilon)\longrightarrow {1\over4},\qquad
 \lambda_1(\mathsf T\mu_\varepsilon)\longrightarrow {4\over9}.
 \tag{2.3}
\]

Indeed, Gaussian convolution gives
\(C_P(X_k+\sqrt\varepsilon G)\leq C_P(X_k)+\varepsilon\), which after
variance normalization gives the lower limits for the gaps.  The reverse
upper limits follow by inserting compactly supported smooth quasimodes for
the Gamma bottom spectral edges; their variances and energies converge.
No spectral attainment is used.

On the other hand,

\[
 J(\mu_\varepsilon)\longrightarrow\infty,\qquad
 J(\mathsf T\mu_\varepsilon)\longrightarrow\infty.
 \tag{2.4}
\]

Fisher information is lower semicontinuous under this convergence.  The
first limiting square-root density has a nonzero jump at its boundary,
while the second behaves as \(t^{1/2}\); neither is in \(W^{1,2}\) across
the boundary.

For \(\mathcal L_J(\nu)=1-1/J(\nu)\), Fisher monotonicity gives a
nonnegative decrement, but

\[
 \mathcal L_J(\mu_\varepsilon)
 -\mathcal L_J(\mathsf T\mu_\varepsilon)
 ={1\over J(\mathsf T\mu_\varepsilon)}
  -{1\over J(\mu_\varepsilon)}\longrightarrow0.
 \tag{2.5}
\]

This disproves a universal charge of \(\log(b/a)\) by the decrement of
\(\mathcal L_J\).  The same example defeats every continuous bounded
transform of \(J\) having a finite limit at infinity.  It is
one-dimensional, so an operator version cannot repair it.

## 3. High Hermite modes destroy fixed-noise MMSE charging

Let \(H_m\) be the probabilists' Hermite polynomial,
\(e_m=H_m/\sqrt{m!}\), and fix even \(m\geq4\).  Put

\[
 V_{\varepsilon,m}(x)={x^2\over2}+\varepsilon e_m(x).
 \tag{3.1}
\]

### 3.1 A nonempty convexity range

Since

\[
 e_m''(x)=\sqrt{m(m-1)}\,e_{m-2}(x),
\]

and \(m-2\) is even, \(e_m''\) has positive leading coefficient and tends
to \(+\infty\) at both ends.  Define

\[
 K_m=\max\{0,-\min_x e_m''(x)\}<\infty.
\]

If \(K_m>0\), every \(0<\varepsilon\leq(2K_m)^{-1}\) satisfies

\[
 V_{\varepsilon,m}''(x)=1+\varepsilon e_m''(x)\geq {1\over2}.
 \tag{3.2}
\]

If \(K_m=0\), (3.2) holds for every positive \(\varepsilon\).  Thus the
admissible interval is nonempty for every fixed \(m\), and strong
convexity proves integrability.

Let \(\widetilde\mu_{\varepsilon,m}\) have density proportional to
\(e^{-V_{\varepsilon,m}}\), and variance-normalize it to
\(\mu_{\varepsilon,m}\).  Symmetry gives zero mean.  Hermite orthogonality
gives

\[
 \operatorname {Var}(\widetilde\mu_{\varepsilon,m})
 =1+m\varepsilon^2+O_m(\varepsilon^3),
\]

so rescaling changes no first-order coefficient:

\[
 {d\mu_{\varepsilon,m}\over d\gamma}
 =1-\varepsilon e_m+O_m(\varepsilon^2).
 \tag{3.3}
\]

### 3.2 Spectral perturbation

In normalized Hermites, the first variations of the mass and energy forms
are

\[
 B_1(f,g)=-\mathbb E_\gamma[fge_m],\qquad
 E_1(f,g)=-\mathbb E_\gamma[f'g'e_m].
\]

At the first Gaussian eigenfunction \(e_1=x\), the only nonzero coupling
after subtracting eigenvalue times mass is

\[
 E_1(e_{m-1},e_1)-B_1(e_{m-1},e_1)=\sqrt m.
 \tag{3.4}
\]

The apparent \(e_{m+1}\) coupling cancels.  Isotropy makes
\(E_\varepsilon(e_1,e_1)=B_\varepsilon(e_1,e_1)=1\) identically, so there
is no diagonal second-order term.  The denominator
\(1-(m-1)=2-m\) therefore gives

\[
 a_{\varepsilon,m}
 =1-{m\over m-2}\varepsilon^2+O_m(\varepsilon^3).
 \tag{3.5}
\]

This one-sided expansion can equivalently be proved by minimizing over
\(e_1+c\varepsilon e_{m-1}\) and controlling the orthogonal complement by
the Gaussian spectral separation.

Normalized self-convolution multiplies the first-order \(m\)-th Hermite
coefficient by

\[
 t_m=2^{1-m/2}.
 \tag{3.6}
\]

The Hermite addition formula proves (3.6).  The same form calculation gives

\[
 b_{\varepsilon,m}
 =1-{m\over m-2}t_m^2\varepsilon^2+O_m(\varepsilon^3).
 \tag{3.7}
\]

Second-order density changes do not alter the displayed coefficient:
isotropy makes the diagonal mass and energy of \(e_1\) equal, while the
off-diagonal correction uses only the first variation.

### 3.3 Gaussian-channel attenuation

For an isotropic \(\nu\), define

\[
 \mathcal F(\nu)=I\left({Z+G\over\sqrt2}\middle\Vert\gamma\right)
 =1-2\,\operatorname {mmse}(Z\mid Z+G),\qquad Z\sim\nu.
 \tag{3.8}
\]

The density ratio in (3.8) is the Ornstein-Uhlenbeck image at correlation
\(2^{-1/2}\).  Thus (3.3) is attenuated to
\(1-\varepsilon2^{-m/2}e_m+O_m(\varepsilon^2)\).  Since
\(\int(e_m')^2d\gamma=m\),

\[
 \mathcal F(\mu_{\varepsilon,m})
 =m2^{-m}\varepsilon^2+O_m(\varepsilon^3),
\]

\[
 \mathcal F(\mathsf T\mu_{\varepsilon,m})
 =m2^{-m}t_m^2\varepsilon^2+O_m(\varepsilon^3).
 \tag{3.9}
\]

Combining (3.5), (3.7), and (3.9) gives

\[
 \boxed{
 \lim_{\varepsilon\downarrow0}
 {\log(b_{\varepsilon,m}/a_{\varepsilon,m})
  \over
  \mathcal F(\mu_{\varepsilon,m})
  -\mathcal F(\mathsf T\mu_{\varepsilon,m})}
 ={2^m\over m-2}.}
 \tag{3.10}
\]

Taking even \(m\to\infty\) disproves a universal fixed-noise
Fisher/MMSE charging inequality inside smooth strongly log-concave
one-dimensional laws.

For comparison,

\[
 \delta(\mu_{\varepsilon,m})={\varepsilon^2\over2}+O_m(\varepsilon^3),
 \qquad
 \delta(\mathsf T\mu_{\varepsilon,m})
 ={t_m^2\varepsilon^2\over2}+O_m(\varepsilon^3),
\]

so the entropy ratio tends \(2m/(m-2)\).  Entropy does not suffer the
fixed-noise high-chaos attenuation.

## 4. Scalar entropy and dependence have no universal product budget

Entropy deficit and mutual information are additive under products.

* For \(\operatorname {Exp}_{\rm iso}^{\otimes n}\), the gap change is
  still \(1/4\mapsto4/9\), while entropy deficits and Hadamard mutual
  informations are \(n\) times their one-dimensional values.  The raw
  budget is not universal, and bounded transforms saturate as
  \(n\to\infty\).
* For
  \(\operatorname {Exp}_{\rm iso}\otimes\gamma^{\otimes(n-1)}\), the same
  gap change occurs but all entropy/dependence production is confined to
  one coordinate.  Division by \(n\) makes the proposed charge tend to
  zero.

Thus neither trace normalization nor the unnormalized scalar quantity is
compatible with both tensorization and Gaussian padding.  A genuinely
directional dependence quantity would be needed.

## 5. Gaussian-channel mean gradients

Let \(X\sim\mu\) be centered, isotropic, and log-concave in
\(\mathbb R^n\), let \(G\sim N(0,I)\) be independent, and put

\[
 Y=X+G,\qquad q=\mathcal L(Y),\qquad
 h(y)=\mathbb E[g(X)\mid Y=y].
\]

Assume \(g\perp\{1,X_1,\ldots,X_n\}\).

### 5.1 The output term

Let \(\rho=\nabla\log q\) and \(J(q)=\mathbb E[\rho\rho^T]\).  The score
identity \(\rho(Y)=-\mathbb E[G\mid Y]\) gives \(J(q)\preceq I\).
Since \(\operatorname {Cov}(Y)=2I\), Cramer-Rao gives
\(J(q)\succeq\frac12I\).  For \(r(Y)=\rho(Y)+Y/2\),

\[
 \mathbb E[rr^T]=J(q)-{1\over2}I\preceq {1\over2}I.
 \tag{5.1}
\]

Affine orthogonality passes through the channel:

\[
 \mathbb Eh=0,\qquad \mathbb E[Yh(Y)]=\mathbb E[Xg(X)]=0.
\]

Integration by parts and conditional Jensen therefore give

\[
 \boxed{
 |\mathbb E_q\nabla h|
 \leq {1\over\sqrt2}\|h\|_{L^2(q)}
 \leq {1\over\sqrt2}\|g\|_{L^2(\mu)}.}
 \tag{5.2}
\]

### 5.2 Exact posterior decomposition

The posterior \(\mu_y=\mathcal L(X\mid Y=y)\) is \(1\)-strongly
log-concave.  Put

\[
 m_y=\mathbb E_yX,\quad A_y=\operatorname {Cov}_y(X),\quad
 a_y=\mathbb E_y\nabla g,
\]

\[
 r_y(x)=g(x)-\mathbb E_yg-a_y\cdot(x-m_y).
 \tag{5.3}
\]

Posterior Poincare applied componentwise to \(\nabla g\) gives

\[
 \mathbb E_y|\nabla r_y|^2
 =\operatorname {Var}_y(\nabla g)
 \leq\mathbb E_y\|D^2g\|_{\rm HS}^2.
 \tag{5.4}
\]

A second posterior Poincare inequality, \(A_y\preceq I\), and
Cauchy-Schwarz give

\[
 |\operatorname {Cov}_y(r_y,X)|
 \leq\big(\mathbb E_y\|D^2g\|_{\rm HS}^2\big)^{1/2}.
 \tag{5.5}
\]

Exponential-family differentiation gives

\[
 \nabla h(y)=\operatorname {Cov}_y(g,X)
 =A_ya_y+\operatorname {Cov}_y(r_y,X).
 \tag{5.6}
\]

Since \(\mathbb Ea_Y=\mathbb E_\mu\nabla g\),

\[
 \boxed{
 \mathbb E_\mu\nabla g-\mathbb E_q\nabla h
 =\mathbb E[(I-A_Y)a_Y]
  -\mathbb E\operatorname {Cov}_Y(r_Y,X).}
 \tag{5.7}
\]

The second term satisfies

\[
 \left|\mathbb E\operatorname {Cov}_Y(r_Y,X)\right|
 \leq\|D^2g\|_{L^2(\mu)}.
 \tag{5.8}
\]

Thus the proposed Hessian estimate is equivalent, up to the additive
constant one in (5.8), to

\[
 \boxed{
 \left|\mathbb E[(I-A_Y)a_Y]\right|
 \leq C\|D^2g\|_{L^2(\mu)}.}
 \tag{AFF}
\]

### 5.3 The affine-mode obstruction

An affine function inside one posterior has zero Hessian but contributes
\((I-A_y)a_y\), generally nonzero.  This explicitly disproves any
pointwise posterior version of (AFF); global affine orthogonality cannot
be imposed posterior by posterior.

The average coefficient is not small:

\[
 \mathbb EA_Y=I-J(q),\qquad
 \mathbb E(I-A_Y)=J(q)\succeq {1\over2}I.
 \tag{5.9}
\]

There is a further exact derivative estimate.  Since
\(Da(y)=\operatorname {Cov}_y(\nabla g,X)\), posterior Poincare and
\(A_y\preceq I\), applied row by row, yield

\[
 \int\|Da(y)\|_{\rm HS}^2q(y)dy
 \leq\int\|D^2g(x)\|_{\rm HS}^2d\mu(x).
 \tag{5.10}
\]

Turning (5.10) into (AFF) requires a global coercive inequality for the
output \(q\).  The known bound is only
\(C_P(q)\leq C_P(\mu)+1\); inserting a dimension-free output Poincare
bound is circular.

KLS would imply the proposed estimate directly.  Two Poincare
applications to affine-orthogonal \(g\) give

\[
 |\mathbb E\nabla g|+\|g\|_2
 \leq C_P(\mu)\|D^2g\|_2.
\]

Together with (5.2), this bounds the left side of (5.7).  Conversely, the
proposed Hessian estimate plus (5.2) gives the desired mean-gradient
graph-norm estimate.  Hence (AFF) is a live KLS-strength gate, not a
consequence of posterior strong convexity or the output Fisher bound.

No hard-support counterexample to the globally orthogonal statement was
found.  Intervals and cubes reduce to the established product estimate,
and the isotropic-ball near-linear residual has left side \(O(n^{-1})\)
versus Hessian norm \(O(n^{-1/2})\).  The audit rules out the proposed
conditional proof, not the global statement (AFF).

## 6. Reverse conditional projection fails on an interval

There is a separate obstruction to reversing the self-convolution argument
by projecting an old slow mode onto the normalized-sum coordinate.  Let
\(f\) be a centered, normalized exact \(\mu\)-eigenfunction with eigenvalue
\(a\), and, for independent \(X,Y\sim\mu\), put

\[
 S={X+Y\over\sqrt2},\qquad D={X-Y\over\sqrt2},
\]

\[
 H(S,D)={f(X)+f(Y)\over\sqrt2},\qquad
 g(S)=\mathbb E[H\mid S],\qquad r=H-g.
 \tag{6.1}
\]

Conditional expectation makes \(g\perp r\), while the product eigenfunction
identity gives

\[
 \mathcal E_{\mu^2}(H,r)=a\langle H,r\rangle
 =a\|r\|_2^2.
\]

Since \(\|H\|_2=1\) and \(\mathcal E_{\mu^2}(H)=a\), expansion of the last
two identities yields the exact formula

\[
 \boxed{\quad
 \mathcal E_{\mathsf T\mu}(g)
 =a(1-2\|r\|_2^2)+\mathcal E_{\mu^2}(r).
 \quad}                                               \tag{6.2}
\]

In particular, there is no general \(H^1\)-contraction under this
conditional expectation.  This already fails for
\(\mu=\operatorname {Unif}[-1,1]\).  Its normalized first nonconstant
Neumann eigenfunction and gap are

\[
 f(x)=\sqrt2\sin {\pi x\over2},\qquad a={\pi^2\over4}.
 \tag{6.3}
\]

Write \(t=\sqrt2s\).  Given \(X+Y=t\), the variable \(X\) is uniform on
an interval of length \(2-|t|\).  Averaging (6.1) on that interval gives

\[
 g(s)={4\over\pi}{\sin(\pi t/2)\over2-|t|},
 \qquad -2<t<2,                                      \tag{6.4}
\]

with the continuous endpoint values.  Define

\[
 J:=\int_0^\pi {\sin^2v\over v}\,dv.
\]

The density of \(t=X+Y\) is \((2-|t|)/4\).  Direct substitution in
(6.4), followed by \(v=\pi(2-t)/2\) on the positive half, gives

\[
 \operatorname {Var}_\nu(g)={8J\over\pi^2},\qquad
 \mathcal E_\nu(g)
 =4\int_0^\pi{(\sin v-v\cos v)^2\over v^3}\,dv
 =4\left(J-{1\over2}\right).                        \tag{6.5}
\]

For the last equality, set \(w(v)=\sin(v)/v\), use
\(w''+2w'/v+w=0\), and integrate
\(\int_0^\pi v(w')^2\,dv\) by parts.  The boundary term vanishes and
\(\int_0^\pi ww'=-1/2\).

The strict comparison in (6.5) has a fully elementary proof.  Split the
integral defining \(J\) at \(\pi/2\).  On the first half,
\(\sin v\geq v-v^3/6\), and hence

\[
 \int_0^{\pi/2}{\sin^2v\over v}\,dv
 \geq {\pi^2\over8}-{\pi^4\over192}
       +{\pi^6\over13824}.                           \tag{6.6}
\]

On the second half put \(u=\pi-v\) and use
\((\pi-u)^{-1}\geq\pi^{-1}+u\pi^{-2}\).  Since

\[
 \int_0^{\pi/2}\sin^2u\,du={\pi\over4},\qquad
 \int_0^{\pi/2}u\sin^2u\,du={\pi^2\over16}+{1\over4},
\]

one obtains

\[
 \int_{\pi/2}^{\pi}{\sin^2v\over v}\,dv
 \geq {5\over16}+{1\over4\pi^2}.                    \tag{6.7}
\]

Consequently

\[
 J\geq {\pi^2\over8}-{\pi^4\over192}
       +{\pi^6\over13824}+{5\over16}+{1\over4\pi^2}
 >{1\over2}+{\pi^2\over16}.                         \tag{6.8}
\]

For completeness, the final strict inequality is also elementary.  If
\(x=\pi^2\), multiplying its difference by \(13824x\) gives

\[
 p(x)=x^4-72x^3+864x^2-2592x+3456.
\]

The polynomial is decreasing on \([9,10]\) and
\(p(x)\geq p(10)=1936>0\); use \(3<\pi<22/7\) to place
\(x\) in that interval.  Equations (6.5) and (6.8) now prove

\[
 \boxed{\quad \mathcal E_\nu(g)>{\pi^2\over4}
 =\mathcal E_{\mu^2}(H).\quad}                       \tag{6.9}
\]

Equivalently, using
\(J=\{\gamma_{\rm E}+\log(2\pi)-\operatorname {Ci}(2\pi)\}/2\),

\[
 \mathcal E_\nu(g)
 =2\{\gamma_{\rm E}+\log(2\pi)-\operatorname {Ci}(2\pi)-1\}
 =2.8753067861\ldots>{\pi^2\over4}.
\]

Thus a reverse argument cannot discard the even-fiber residual in (6.2)
by asserting that conditional expectation contracts Dirichlet energy.

## 7. Hessian refinement of the forward defect estimate

The exact near-equality structure of the valid forward argument can be
made sharper.  Retain the notation of Section 6, but now let \(f\) be a
centered, normalized exact \(\nu=\mathsf T\mu\)-eigenfunction with
eigenvalue \(b\).  Put

\[
 Z=\nabla f(S),\quad U=\mathbb E[Z\mid X],\quad
 V=\mathbb E[Z\mid Y],\quad m=\mathbb EZ,
\]

\[
 W=Z-U-V+m,\qquad v=\operatorname {Var}(U),\qquad
 w=\|W\|_2^2.
\]

For the scalar Hoeffding residual
\(R=f(S)-\mathbb E[f(S)\mid X]-\mathbb E[f(S)\mid Y]\), let
\(E_R=\mathcal E_{\mu^2}(R)\).  Orthogonality of the vector Hoeffding
decomposition gives

\[
 \boxed{\quad
 \operatorname {Var}(Z)=2v+w,\qquad E_R=v+w.
 \quad}                                               \tag{7.1}
\]

There is an elementary lower bound which loses no dimension factor.  The
weak eigenfunction identity tested against the coordinate functions gives

\[
 \mathbb E\nabla f=b\,\mathbb E[Sf].
\]

Because \(\mathbb E|Z|^2=b\) and isotropy implies
\(|\mathbb E[Sf]|\leq1\),

\[
 \boxed{\quad
 \operatorname {Var}(\nabla f)
 =b-b^2|\mathbb E[Sf]|^2\geq b(1-b).
 \quad}                                               \tag{7.2}
\]

Combining (7.1), (7.2), and the scalar defect estimate
\(E_R\leq2(b-a)\) recovers

\[
 b(1-b)\leq\operatorname {Var}(Z)
 \leq2E_R\leq4(b-a).                                 \tag{7.3}
\]

The Hessian identities reveal what saturation of (7.3) would require.
Set

\[
 H=D^2f(S),\quad B=\mathbb E\|H\|_{\rm HS}^2,
 \quad K(X)=\mathbb E[H\mid X],\quad C=\mathbb EH,
\]

\[
 Q=H-K(X)-K(Y)+C.
\]

In the smooth convex setting, the integrated Bochner identity and
\(A_\nu f=bf\) give

\[
 B\leq\|A_\nu f\|_2^2=b^2.                          \tag{7.4}
\]

The matrix-valued Hoeffding decomposition is orthogonal, so

\[
 B=2\mathbb E\|K-C\|_{\rm HS}^2
   +\mathbb E\|Q\|_{\rm HS}^2+\|C\|_{\rm HS}^2.     \tag{7.5}
\]

Moreover, \(DU=K/\sqrt2\).  Applying the \(\mu\)-Poincare inequality
componentwise to \(U\), and then using (7.5), gives

\[
 2av\leq\mathbb E\|K\|_{\rm HS}^2
 \leq {B+\|C\|_{\rm HS}^2\over2},
\]

or

\[
 \boxed{\quad4av\leq B+\|C\|_{\rm HS}^2.\quad}      \tag{7.6}
\]

Finally, condition on \(X\) and apply Poincare in the \(Y\)-variable to
\(Z=\nabla f((X+Y)/\sqrt2)\).  Since
\(D_YZ=H/\sqrt2\),

\[
 \boxed{\quad B\geq2aE_R.\quad}                     \tag{7.7}
\]

Consider now a sequence of such exact modes for which

\[
 a\longrightarrow0,\qquad {b\over a}\longrightarrow{4\over3}.
 \tag{7.8}
\]

The nonlinear extremal regime may additionally impose a vanishing affine
projection \(|\mathbb E[Sf]|\to0\), but the following conclusion only
uses the universal bound \(|\mathbb E[Sf]|\leq1\).  Equality squeezing in
(7.3) and (7.1) yields

\[
 w=o(a),\qquad E_R=v={2a\over3}+o(a).                 \tag{7.9}
\]

Equations (7.4) and (7.8) give

\[
 B\leq{16a^2\over9}+o(a^2),                          \tag{7.10}
\]

whereas (7.6) and (7.9) give

\[
 B+\|C\|_{\rm HS}^2\geq{8a^2\over3}+o(a^2).
 \tag{7.11}
\]

Since (7.5) also gives \(\|C\|_{\rm HS}^2\leq B\), the ratio is
well-defined along every sufficiently far nontrivial mode, and

\[
 \boxed{\quad
 \liminf {\|\mathbb E D^2f(S)\|_{\rm HS}^2
                 \over
                 \mathbb E\|D^2f(S)\|_{\rm HS}^2}
 \geq {1\over2}.
 \quad}                                               \tag{7.12}
\]

Thus a sequence saturating the \(4/3\) forward amplification must place at
least half of its Hessian energy in the mean Hessian.  Excluding precisely
this coherent mean-Hessian, or directional third-moment, concentration is
the remaining rigidity gate; estimates of the oscillatory Hessian alone
cannot do so.
