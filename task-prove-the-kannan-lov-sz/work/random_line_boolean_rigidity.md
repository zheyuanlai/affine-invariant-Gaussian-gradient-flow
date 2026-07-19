# Boolean rigidity for the variance-normalized random-line form

## 0. Verdict

Let \(\mu\) be an isotropic log-concave probability on
\(\mathbb R^n\), let \(E\) be balanced, and put

\[
 h=2\mathbf 1_E-1,
 \qquad
 \mathcal Q_\mu(h)=
 \mathbb E_\theta\int_{\theta^\perp}
 {\operatorname {Var}(h\mid \pi_{\theta^\perp}X=y)
  \over \sigma_{\theta,y}}\,d(\pi_{\theta^\perp}\mu)(y).
 \tag{0.1}
\]

The proposed Boolean-affine rigidity estimate is

\[
 \operatorname {dist}_{L^2(\mu)}(h,\operatorname {Affine})^2
 \le C\sqrt n\,\mathcal Q_\mu(h).                         \tag{0.2}
\]

The main conclusion of this report is an exact obstruction, rather than
a new proof of random-line uncertainty.

**Equivalence theorem.** There is a universal \(c_0>0\) such that every
balanced Boolean \(h\) under every isotropic log-concave \(\mu\) satisfies

\[
 c_0\le
 \operatorname {dist}_{L^2(\mu)}(h,\operatorname {Affine})^2
 \le1.                                                     \tag{0.3}
\]

The self-contained proof below gives \(c_0=1/96\); the sharp
one-dimensional mode--variance bound improves this to \(1/48\).
Consequently, uniformly over this class,
(0.2) holds with a universal constant if and only if

\[
                         \mathcal Q_\mu(h)\ge {c\over\sqrt n}.\tag{0.4}
\]

Moreover, if \(\mathcal U_\mu(E)\) is the random-line uncertainty from
`random_line_uncertainty.md`, then

\[
       2\mathcal U_\mu(E)\le\mathcal Q_\mu(h)
       \le4\mathcal U_\mu(E).                              \tag{0.5}
\]

Thus Boolean-affine rigidity is quantitatively equivalent to the
missing uncertainty theorem.  Removing the affine component does not
weaken the KLS-strength step: balanced Boolean functions have a fixed
amount of nonlinear \(L^2\) mass already.

The restriction to Boolean functions is essential.  For the standard
Gaussian and

\[
 f(x)={|x|^2-n\over\sqrt{2n}},                             \tag{0.6}
\]

one has \(f\perp\operatorname {Affine}\), \(\|f\|_2=1\), but

\[
                         \mathcal Q_{\gamma_n}(f)={1\over n}.\tag{0.7}
\]

Hence (0.2) is false for general nonlinear functions by a factor
\(\sqrt n\).  Taking the sign of (0.6) creates a codimension-one jump;
its line energy becomes \(\Theta(n^{-1/2})\).  Exact Gaussian Hermite
and Laguerre calculations below show that this gain is supplied by an
entire hierarchy of polynomial degrees extending to order \(n\), not by
a second-chaos estimate.

All requested models pass (0.2), at the sharp scale when the interface
is a single cut:

* the Gaussian halfspace has affine residual \(1-2/\pi\) and
  \(\mathcal Q\asymp n^{-1/2}\);
* the cube coordinate cut has affine residual \(1/4\) and
  \(\mathcal Q\asymp n^{-1/2}\), whereas the full parity checkerboard
  has residual \(1\) and \(\mathcal Q\asymp\sqrt n\);
* the Gaussian median sphere has residual \(1\) and
  \(\mathcal Q\asymp n^{-1/2}\);
* a balanced regular-simplex cap has an explicit residual tending to
  \(1-(\log2)^2\) and \(\mathcal Q\asymp n^{-1/2}\);
* for the product-exponential median maximum the affine projection tends
  to zero and \(\mathcal Q\asymp n^{-1/2}\), while parity again has
  \(\mathcal Q\asymp\sqrt n\).

No example in this list disproves (0.2).  The theorem above instead
shows why proving it is not an intermediate lemma: it is the random-line
conjecture with a bounded, automatically nonzero numerator.

## 1. The line form and its Boolean normalization

For a fixed \(\theta\), disintegrate \(\mu\) over
\(\pi_{\theta^\perp}\).  On the fiber over \(y\), write

\[
 p_{\theta,y}=\mathbb P(X\in E\mid\pi_{\theta^\perp}X=y),
 \qquad q_{\theta,y}=\min(p_{\theta,y},1-p_{\theta,y}).    \tag{1.1}
\]

Since \(h\in\{-1,1\}\),

\[
 \operatorname {Var}(h\mid y)=4p_{\theta,y}(1-p_{\theta,y}).\tag{1.2}
\]

For \(0\le p\le1\),

\[
                         2\min(p,1-p)
 \le4p(1-p)\le4\min(p,1-p).                               \tag{1.3}
\]

After division by \(\sigma_{\theta,y}\), integration in \(y\), and
averaging in \(\theta\), this proves (0.5).

It is useful to polarize (0.1).  Define

\[
 \mathcal E_\mu(f,g)=
 \mathbb E_\theta\int {\operatorname {Cov}(f,g\mid y)
                   \over\sigma_{\theta,y}}\,d\mu_\theta(y).
 \tag{1.4}
\]

This is a positive semidefinite symmetric form and
\(\mathcal Q_\mu(f)=\mathcal E_\mu(f,f)\).  For an affine linear
function \(\ell_a(x)=a\cdot x\),

\[
 \mathcal Q_\mu(\ell_a)=
 \mathbb E_\theta (a\cdot\theta)^2
       \int\sigma_{\theta,y}\,d\mu_\theta(y).             \tag{1.5}
\]

Conditional variance decomposition and isotropy give

\[
 \int\sigma_{\theta,y}^2d\mu_\theta(y)\le1,
 \qquad
 \mathcal Q_\mu(\ell_a)\le {|a|^2\over n}.               \tag{1.6}
\]

Thus the affine component naturally lives at scale \(n^{-1}\), below
the desired Boolean scale \(n^{-1/2}\).  If \(h=\ell_a+r\) is the
orthogonal affine decomposition, Cauchy--Schwarz for (1.4) gives only

\[
 \sqrt{\mathcal Q_\mu(r)}
 \le\sqrt{\mathcal Q_\mu(h)}+{|a|\over\sqrt n}.            \tag{1.7}
\]

A hypothetical bad Boolean function therefore leaves a genuinely
nonlinear residual with small line energy.  There is no separately
controlled affine term whose subtraction closes the estimate.

## 2. Universal slack of the best affine predictor

We prove (0.3) with an explicit constant.

### Lemma 2.1 (bounded-density variance)

If a probability density \(g\) on \(\mathbb R\) satisfies
\(\|g\|_\infty\le L\), then

\[
                         \operatorname {Var}(g)
 \ge {1\over12L^2}.                                       \tag{2.1}
\]

#### Proof

For any \(b\in\mathbb R\), the bathtub principle says that among
densities bounded by \(L\) and having total mass one, the integral
\(\int(x-b)^2g(x)dx\) is minimized by putting density \(L\) on the
interval of length \(1/L\) centered at \(b\).  The minimum is
\(1/(12L^2)\).  Minimize the left side over \(b\).  This proves (2.1).
\(\square\)

We next record a self-contained, slightly nonsharp density bound.  If
\(f\) is a log-concave probability density with mode \(x_0\), height
\(M=f(x_0)\), and variance \(\sigma^2\), then

\[
                         M\le {\sqrt2\over\sigma}.          \tag{2.2}
\]

Indeed, after translating \(x_0=0\), put
\(A=\int_0^\infty f\) and \(B=1-A\).  The survival function
\(S_+(x)=\int_x^\infty f\) is log-concave by Prekopa--Leindler.  Since
\(S_+(0)=A\) and its right logarithmic derivative at zero is
\(-M/A\), concavity gives

\[
                         S_+(x)\le A e^{-Mx/A}\quad(x\ge0).\tag{2.3}
\]

The analogous bound holds on the negative half-line with mass \(B\).
Integration by parts yields

\[
 \int_0^\infty x^2f(x)dx\le {2A^3\over M^2},
 \qquad
 \int_{-\infty}^0 x^2f(x)dx\le {2B^3\over M^2}.           \tag{2.4}
\]

Therefore

\[
 \sigma^2\le\mathbb E(X-x_0)^2
 \le {2(A^3+B^3)\over M^2}\le {2\over M^2},
\]

which proves (2.2).  The sharp one-dimensional estimate is

\[
 \|f\|_\infty\le {1\over\sqrt{\operatorname {Var}(\nu)}}.\tag{2.5}
\]

This is the upper half of the one-dimensional mode--variance estimate
\(1/(12\operatorname {Var}\nu)\le\|f\|_\infty^2
\le1/\operatorname {Var}\nu\).  We do not need its sharp constant: the
proved estimate (2.2) is enough below.

### Lemma 2.2 (affine slack)

Let \(\mu\) be isotropic and log-concave, and let
\(h:\mathbb R^n\to\{-1,1\}\) satisfy \(\int h\,d\mu=0\).  Then

\[
 \inf_{b\in\mathbb R,\,a\in\mathbb R^n}
 \int(h-b-a\cdot x)^2d\mu(x)
 \ge {1\over96}.                                          \tag{2.6}
\]

#### Proof

Since \(1,x_1,\ldots,x_n\) are orthogonal in \(L^2(\mu)\), the best
affine predictor is \(a\cdot x\), where

\[
                         a=\int xh(x)d\mu(x),              \tag{2.7}
\]

and the squared residual is \(1-|a|^2\).  If \(a=0\), there is nothing
to prove.  Otherwise put \(u=a/|a|\) and \(Y=u\cdot X\).  The law of
\(Y\) is a one-dimensional isotropic log-concave probability with a
continuous density.  Let \(m\) be its median and put
\(s(Y)=\operatorname {sign}(Y-m)\).  If
\(g(Y)=\mathbb E(h\mid Y)\), then \(|g|\le1\) and \(\mathbb Eg=0\).
Pointwise,

\[
 (Y-m)(s(Y)-g(Y))\ge0.                                    \tag{2.8}
\]

Since both \(s\) and \(g\) have mean zero, (2.5) implies

\[
 |a|=\mathbb E(Yh)=\mathbb E(Yg(Y))
 \le\mathbb E(Ys(Y))=\mathbb E|Y-m|.                      \tag{2.9}
\]

Let \(v_+\) and \(v_-\) be the conditional variances of \(Y\) on the
upper and lower median half-lines.  By (2.2), each conditional density
is bounded by \(2\|f_Y\|_\infty\le2\sqrt2\).  Lemma 2.1 gives

\[
                         v_+,v_-\ge {1\over96}.             \tag{2.10}
\]

If \(m_+\) and \(m_-\) are the two conditional means, then
\(m_-=-m_+\), because the halves have mass \(1/2\) and \(\mathbb EY=0\).
Furthermore,

\[
 \mathbb E|Y-m|=\mathbb E(Ys(Y))=m_+.
\]

The law of total variance therefore gives

\[
 1=m_+^2+{v_++v_-\over2},
 \qquad
 1-|a|^2\ge1-m_+^2={v_++v_-\over2}\ge{1\over96}.          \tag{2.11}
\]

This is (2.3). \(\square\)

### Theorem 2.3 (rigidity is uncertainty)

For balanced Boolean functions under isotropic log-concave measures,
the following statements are quantitatively equivalent:

1. there is a universal \(c>0\) such that
   \(\mathcal Q_\mu(h)\ge c/\sqrt n\);
2. there is a universal \(C<\infty\) such that
   \(\operatorname {dist}(h,\operatorname {Affine})^2
      \le C\sqrt n\,\mathcal Q_\mu(h)\).

#### Proof

Statement 1 implies statement 2 because the squared distance is at most
\(\|h\|_2^2=1\).  Statement 2 and Lemma 2.2 imply

\[
 \mathcal Q_\mu(h)\ge {1\over96C\sqrt n}.
\]

\(\square\)

Combining Theorem 2.3 with (0.5) proves the claimed equivalence with the
random-line uncertainty estimate.  In particular, a proof of (0.2)
would complete the same KLS-strength step as a direct proof of
\(\mathcal U_\mu(E)\ge c/\sqrt n\).

## 3. Exact Gaussian operator calculations

Let \(\gamma_n\) be standard Gaussian measure.  Here every conditional
line variance is one.  Let

\[
 P_\theta f=\mathbb E(f(X)\mid\pi_{\theta^\perp}X),
 \qquad K_n=\mathbb E_\theta P_\theta.                     \tag{3.1}
\]

Since \(P_\theta\) is an orthogonal projection,

\[
 \mathcal Q_{\gamma_n}(f)
 =\langle f,(I-K_n)f\rangle_{L^2(\gamma_n)}.               \tag{3.2}
\]

The operator \(K_n\) preserves every Wiener chaos and commutes with the
orthogonal group.  The following exact multipliers exhibit the slowest
trace-heavy modes.

### 3.1 First and second chaos

For \(f(x)=a\cdot x\),

\[
                         K_nf=\left(1-{1\over n}\right)f.  \tag{3.3}
\]

For the unit radial quadratic

\[
 R_1(x)={|x|^2-n\over\sqrt{2n}},
\]

conditional expectation deletes the centered square in the sampled
direction.  Hence

\[
 K_nR_1={n-1\over n}R_1,
 \qquad \mathcal Q_{\gamma_n}(R_1)={1\over n}.             \tag{3.4}
\]

If \(A=A^T\), \(\operatorname {tr}A=0\), and
\(T_A(x)=x^TAx\), then averaging

\[
 \|P_{\theta^\perp}AP_{\theta^\perp}\|_F^2
 =\|A\|_F^2-2|A\theta|^2+(\theta^TA\theta)^2
\]

gives the traceless-quadratic eigenvalue

\[
 {\langle T_A,K_nT_A\rangle\over\|T_A\|_2^2}
 =1-{2(n+1)\over n(n+2)}.                                 \tag{3.5}
\]

Thus removing affine functions does not create an
\(n^{-1/2}\)-spectral gap: even the radial second chaos has gap exactly
\(1/n\).

### 3.2 Radial even chaoses

Put \(a=n/2\), \(S=|X|^2/2\), and let

\[
 R_j(x)=c_{n,j}L_j^{a-1}(S),                               \tag{3.6}
\]

where \(L_j^\alpha\) is the generalized Laguerre polynomial and
\(c_{n,j}\) makes \(R_j\) a unit vector.  This spans the radial part of
chaos \(2j\).  The Laguerre generating function and
\(T^2/2\sim\operatorname {Gamma}(1/2,1)\) give

\[
 \mathbb E_TL_j^{a-1}(S_\perp+T^2/2)
 =L_j^{a-3/2}(S_\perp).                                   \tag{3.7}
\]

Using

\[
 \mathbb E_{G\sim\operatorname {Gamma}(b,1)}
       [L_j^{b-1}(G)^2]={ (b)_j\over j!},                 \tag{3.8}
\]

we obtain the exact eigenvalue

\[
 \boxed{\quad
 K_nR_j=\kappa^{\rm rad}_{n,j}R_j,
 \qquad
 \kappa^{\rm rad}_{n,j}=
 {((n-1)/2)_j\over(n/2)_j}
 =\prod_{r=0}^{j-1}\left(1-{1\over n+2r}\right).
 \quad}                                                    \tag{3.9}
\]

In particular, for \(n\ge2\),

\[
 1-\kappa^{\rm rad}_{n,j}\asymp {j\over n}
 \quad(1\le j\le n),
 \qquad
 \kappa^{\rm rad}_{n,j}\asymp\sqrt{n\over n+j}
 \quad(j\ge n),                                          \tag{3.10}
\]

with universal implicit constants.

### 3.3 The trace-heavy odd chaoses

For a unit vector \(v\), the functions

\[
 V_{j,v}(x)=d_{n,j}(v\cdot x)L_j^{a}(S)                    \tag{3.11}
\]

form the standard-vector irreducible component of chaos \(2j+1\).
Conditioning in direction \(\theta\) kills the odd-in-\(T\) term and
replaces \(v\) by \(P_{\theta^\perp}v\).  Since

\[
 \mathbb E[(v\cdot X)^2L_j^a(S)^2]={ (a+1)_j\over j!},    \tag{3.12}
\]

Schur's lemma and \(\mathbb E|P_{\theta^\perp}v|^2=1-1/n\)
give

\[
 \boxed{\quad
 K_nV_{j,v}=\kappa^{\rm vec}_{n,j}V_{j,v},
 \qquad
 \kappa^{\rm vec}_{n,j}=
 \left(1-{1\over n}\right)
 {((n+1)/2)_j\over((n+2)/2)_j}.
 \quad}                                                    \tag{3.13}
\]

Again \(1-\kappa^{\rm vec}_{n,j}\asymp(j+1)/n\) for
\(j\le n\).  Trace contraction therefore makes polynomial degree
\(2j\) behave like only \(j\) random-direction constraints.  A proof
based on finitely many low chaoses cannot produce the Boolean square-root
gain.

### 3.4 Rank-one Hermite rays

Let \(H_k(v\cdot x)/\sqrt{k!}\) be a unit rank-one vector in chaos
\(k\).  If \(b=v\cdot\theta\), the conditional projection has squared
norm \((1-b^2)^k\).  Since \(b^2\) is
\(\operatorname {Beta}(1/2,(n-1)/2)\), its exact Rayleigh multiplier is

\[
 \kappa^{\rm ray}_{n,k}
 =\mathbb E(1-b^2)^k
 ={((n-1)/2)_k\over(n/2)_k}.                              \tag{3.14}
\]

This is an eigenvalue only when the corresponding rank-one tensor lies
in a single irreducible component; in general (3.14) is the exact
Rayleigh quotient.  The distinction is important, but the formula is
enough for the halfspace calculation below.

## 4. Gaussian halfspace and the Boolean degree hierarchy

Let \(h(x)=\operatorname {sign}(x_1)\).  Its best affine coefficient is

\[
                         a=\sqrt{2\over\pi}\,e_1,
 \qquad
 \operatorname {dist}(h,\operatorname {Affine})^2
 =1-{2\over\pi}.                                          \tag{4.1}
\]

For fixed \(\theta\), put \(b=|\theta_1|\).  Two independent samples
on the conditional line have first-coordinate correlation \(1-b^2\).
Sheppard's formula gives

\[
 \mathcal Q_{\gamma_n,\theta}(h)
 ={2\over\pi}\arccos(1-b^2).                             \tag{4.2}
\]

Consequently

\[
                         \mathcal Q_{\gamma_n}(h)
 ={2\over\pi}\mathbb E_\theta\arccos(1-\theta_1^2)
 \asymp\mathbb E|\theta_1|\asymp {1\over\sqrt n}.        \tag{4.3}
\]

The same conclusion can be read degree by degree.  The normalized
Hermite weights of the sign function are

\[
 w_{2m+1}={2\over\pi}
 {\binom{2m}{m}\over4^m(2m+1)},
 \qquad w_{2m}=0,                                         \tag{4.4}
\]

and \(w_{2m+1}\asymp(m+1)^{-3/2}\).  Equations (3.14) and
(4.4) yield the exact identity

\[
 \mathcal Q_{\gamma_n}(h)
 =\sum_{m\ge0}w_{2m+1}
       \left(1-kappa^{\rm ray}_{n,2m+1}\right).           \tag{4.5}
\]

For \(k\le n\), the multiplier loss is \(\asymp k/n\); for
\(k\ge n\), it is bounded below once \(k/n\) is bounded below.  Hence

\[
 {1\over n}\sum_{k\le n}k w_k
 +\sum_{k>n}w_k\asymp {1\over\sqrt n}.                    \tag{4.6}
\]

The square-root gain is therefore not located in one nonlinear mode.
It is the accumulation of the \(k^{-3/2}\) Boolean jump tail through
degrees \(k\asymp n\).

This calculation also identifies the burden of a polynomial-mode proof
for a general log-concave measure.  A generic spectral estimate after
affine projection can give only a cost of order \(k/n\) to a mode of
effective degree \(k\).  To reach \(n^{-1/2}\), one needs a sharp
Boolean approximation statement at all degrees up to \(n\), with the
one-jump tail as the extremal case.  A fixed-degree or bounded-depth
hierarchy stops at order \(1/n\).

## 5. Gaussian radial quadratic threshold

Let \(s_n\) be the median of \(S=|X|^2/2\sim
\operatorname {Gamma}(n/2,1)\), and put

\[
                         h_{\rm rad}(x)
 =2\mathbf1_{\{S\le s_n\}}-1.                             \tag{5.1}
\]

Rotational symmetry and evenness imply

\[
 \int xh_{\rm rad}(x)d\gamma_n(x)=0,
 \qquad
 \operatorname {dist}(h_{\rm rad},\operatorname {Affine})^2=1.\tag{5.2}
\]

For any \(\theta\), write

\[
 A={|P_{\theta^\perp}X|^2\over2}
 \sim\operatorname {Gamma}((n-1)/2,1),
 \qquad B={\langle X,\theta\rangle^2\over2}
 \sim\operatorname {Gamma}(1/2,1).                        \tag{5.3}
\]

They are independent and the fiber probability is
\(p(A)=\mathbb P(B\le s_n-A)\), interpreted as zero when \(A>s_n\).
Thus

\[
 \mathcal Q_{\gamma_n}(h_{\rm rad})
 =4\mathbb E[p(A)(1-p(A))].                               \tag{5.4}
\]

On the fixed window \(s_n-1\le A\le s_n-1/4\), both \(p(A)\) and
\(1-p(A)\) are bounded below universally.  The gamma density of \(A\)
on that window is \(\asymp n^{-1/2}\), by Stirling's inequalities and
the fact that \(s_n=n/2+O(1)\).  This gives the lower bound in

\[
                         \mathcal Q_{\gamma_n}(h_{\rm rad})
 \asymp {1\over\sqrt n}.                                  \tag{5.5}
\]

For the upper bound, the density of
\(A\sim\operatorname {Gamma}((n-1)/2,1)\) is everywhere at most
\(C/\sqrt n\) when \(n\ge3\).  With \(r=s_n-A\),
\[
 \mathbb E[p(A)(1-p(A))]
 \le {C\over\sqrt n}\int_0^\infty
       F_{1/2}(r)(1-F_{1/2}(r))\,dr.                       \tag{5.5a}
\]
The integral is finite: near zero \(F_{1/2}(r)=O(\sqrt r)\), and at
infinity \(1-F_{1/2}(r)=O(r^{-1/2}e^{-r})\).  The dimensions \(n=1,2\)
are absorbed by changing the universal constants.  This proves the
other direction of (5.5).

Expanding (5.1) in the radial Laguerre basis gives

\[
 h_{\rm rad}=\sum_{j\ge1}b_{n,j}R_j,
 \qquad
 \mathcal Q_{\gamma_n}(h_{\rm rad})
 =\sum_{j\ge1}b_{n,j}^2(1-\kappa^{\rm rad}_{n,j}).         \tag{5.6}
\]

The smooth statistic \(R_1\) alone costs \(1/n\), while thresholding it
creates the full Laguerre tail needed for (5.5).  Thus the radial model
rules out any argument claiming that affine orthogonality plus a
second-order mode estimate is enough.

## 6. Product and convex-body models

Throughout this section, (0.5) transfers the estimates for
\(\mathcal U\) proved in `random_line_uncertainty.md` to
\(\mathcal Q\).

### 6.1 Cube: coordinate cut and Walsh parity

Let \(\mu\) be uniform on \([ -\sqrt3,\sqrt3]^n\) and
\(h(x)=\operatorname {sign}(x_1)\).  Then

\[
 a_1=\mathbb E|X_1|={\sqrt3\over2},
 \qquad a_j=0\ (j>1),
 \qquad
 \operatorname {dist}(h,\operatorname {Affine})^2={1\over4}.\tag{6.1}
\]

The exact chord calculation gives

\[
 {a_n\over48}\le\mathcal U_\mu(E)\le {a_n\over2},
 \qquad a_n=\mathbb E|\theta_1|\asymp n^{-1/2}.            \tag{6.2}
\]

Therefore

\[
                         {a_n\over24}
 \le\mathcal Q_\mu(h)\le2a_n.                            \tag{6.3}
\]

Now let \(\varepsilon_i=\operatorname {sign}(x_i)\) and
\(h_{\rm par}=\prod_{i=1}^n\varepsilon_i\).  For \(n\ge2\), every
constant and affine coefficient vanishes: in Walsh language this is a
pure degree-\(n\) threshold bit.  Hence its affine residual is one.  The
checkerboard calculation gives

\[
                         c\sqrt n\le\mathcal Q_\mu(h_{\rm par})
 \le C\sqrt n.                                            \tag{6.4}
\]

Thus high Walsh oscillation makes the inequality easier; the coordinate
cut, not parity, is the sharp cube test.

### 6.2 Balanced regular-simplex cap

Use barycentric coordinates \(\lambda_0,\ldots,\lambda_n\), normalized
as in `random_line_uncertainty.md`, and let

\[
 E=\{\lambda_0\ge t_n\},
 \qquad t_n=1-2^{-1/n},
 \qquad h=2\mathbf1_E-1.                                  \tag{6.5}
\]

By the symmetry fixing vertex \(v_0\), the affine coefficient is in the
direction \(u=v_0/|v_0|\).  Since

\[
 u\cdot x=\sqrt{{n+2}\over n}\big((n+1)\lambda_0-1\big)  \tag{6.6}
\]

and \(\lambda_0\sim\operatorname {Beta}(1,n)\), direct integration
gives

\[
                         |a|=\sqrt{n(n+2)}\,t_n.            \tag{6.7}
\]

Consequently

\[
 \operatorname {dist}(h,\operatorname {Affine})^2
 =1-n(n+2)(1-2^{-1/n})^2
 \longrightarrow1-(\log2)^2.                             \tag{6.8}
\]

The simplex chord-race estimate gives

\[
                         \mathcal Q_\mu(h)\asymp n^{-1/2}.\tag{6.9}
\]

### 6.3 Product exponentials: maximum and parity

Let \(X_i=Z_i-1\), where the \(Z_i\) are independent unit exponentials.
For the median maximum, put

\[
 E_{\max}=\{\max_i Z_i\le H_n\},
 \quad q_n=e^{-H_n},\quad d_n=1-q_n,\quad d_n^n={1\over2},
 \quad h_{\max}=2\mathbf1_{E_{\max}}-1.                   \tag{6.10}
\]

All affine coefficients are equal.  Since

\[
 \int_0^{H_n}(z-1)e^{-z}dz=-H_nq_n,                       \tag{6.11}
\]

we obtain exactly

\[
 a_i=-{H_nq_n\over d_n},
 \qquad
 \operatorname {dist}(h_{\max},\operatorname {Affine})^2
 =1-n\left({H_nq_n\over d_n}\right)^2.                   \tag{6.12}
\]

Here \(q_n\sim(\log2)/n\) and \(H_n=\log n+O(1)\), so the residual
tends to one.  The facet-race calculation gives

\[
                         \mathcal Q_\mu(h_{\max})
 \asymp n^{-1/2}.                                         \tag{6.13}
\]

For median parity,

\[
 h_{\rm par}=\prod_{i=1}^n
       \operatorname {sign}(Z_i-\log2).                   \tag{6.14}
\]

When \(n\ge2\), independence and the zero mean of every sign bit imply
\(\mathbb E(X_i h_{\rm par})=0\) for all \(i\).  Thus the affine
residual is one, while

\[
                         \mathcal Q_\mu(h_{\rm par})
 \asymp\sqrt n.                                           \tag{6.15}
\]

As for the cube, the many-interface example is far from extremal.

## 7. What a polynomial hierarchy would have to prove

The Gaussian formulas make the required hierarchy quantitative.  Write
schematically

\[
 h=\ell+\sum_{k\ge2}h_k,
 \qquad W_k=\|h_k\|_2^2.                                  \tag{7.1}
\]

Trace-heavy modes show that the strongest possible generic spectral
cost at effective degree \(k\le n\) is only \(k/n\).  A Boolean proof
must therefore supply high-degree mass.  If a tail estimate of the
sharp one-jump form

\[
                         \sum_{k>d}W_k\gtrsim d^{-1/2}      \tag{7.2}
\]

is available through \(d\asymp n\), then the degrees above \(d\) cost
at least \(d/n\), and optimization gives

\[
 {d\over n}\sum_{k>d}W_k
 \gtrsim {\sqrt d\over n},
 \qquad d\asymp n
 \quad\Longrightarrow\quad n^{-1/2}.                     \tag{7.3}
\]

By contrast, a \(d^{-1}\) approximation tail yields only \(n^{-1}\),
and any hierarchy stopped at dimension-free degree also yields only
\(n^{-1}\).  Halfspaces show that the exponent \(1/2\) and depth
\(d\asymp n\) are both sharp.

For a general log-concave measure there is no orthogonal-group
diagonalization of the line operator, and no measure-independent notion
of Hermite degree.  Replacing (7.2) by the intrinsic assertion that a
balanced Boolean function has line energy at least \(n^{-1/2}\) is
exactly (0.4).  Theorem 2.3 shows that an affine residual on the left does
not reduce this burden.

One can phrase this as a clean no-go statement.

### Corollary 7.1 (a precise non-Boolean no-go)

There is no universal \(C\) for which

\[
 \operatorname {dist}_{L^2(\gamma_n)}
      (f,\operatorname {Affine})^2
 \le C\sqrt n\,\mathcal Q_{\gamma_n}(f)                   \tag{7.4}
\]

holds for all \(n\) and all \(f\in L^2(\gamma_n)\).  Indeed, the unit
radial quadratic (0.6) makes the two sides equal to \(1\) and
\(C/\sqrt n\), respectively.  More generally, deleting any fixed
number of low radial chaoses still leaves a unit radial chaos with line
gap \(O(1/n)\), by (3.9).

Thus a valid proof for Boolean functions must use the jump constraint;
a generic higher-order spectral gap after affine projection is false.
The halfspace calculation further shows that the natural sharp Boolean
mechanism uses effective degrees through order \(n\), or an equivalent
non-spectral global mechanism.

This does not disprove Boolean rigidity.  It identifies it as the full
random-line uncertainty statement and rules out the proposed affine
decomposition as an independent shortcut.

## 8. Scope and rigor notes

1. The full-dimensional assumption is harmless here.  A log-concave
   measure supported on a \(k\)-dimensional affine space is identified
   isometrically with that space, and all occurrences of \(n\),
   isotropy, directions, and affine functions are then interpreted in
   dimension \(k\).
2. Projections of a full-dimensional log-concave probability have
   continuous log-concave densities, so the median split in Lemma 2.2
   has no atom.  The same proof on the affine support covers degenerate
   ambient covariance.
3. The only general log-concave input in the equivalence theorem is the
   one-dimensional mode--variance bound (2.2).  No Poincare inequality
   for arbitrary functions is assumed.
4. The Gaussian spectral calculations are model audits, not inputs to
   Theorem 2.3.  In particular, no Gaussian log-Sobolev statement is
   transferred to a general log-concave measure.
