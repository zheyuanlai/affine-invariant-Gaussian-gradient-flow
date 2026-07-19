# Boolean fibers, hit-and-run mixing, and Beckmann transport

## 0. Outcome

Let \(\mu\) be an isotropic log-concave probability on \(\mathbb R^n\),
let \(E\) have mass \(1/2\), and put

\[
 q=1_E-1_{E^c}.
\]

For a unit vector \(\theta\), disintegrate \(\mu\) over the lines parallel
to \(\theta\), write

\[
 p_\theta(y)=\mu(E\mid P_{\theta^\perp}X=y),\qquad
 I_\theta=\mathbb E_y\min(p_\theta(y),1-p_\theta(y)).          \tag{0.1}
\]

There are four exact conclusions.

1.  \(I_\theta\) is, up to a factor two, the probability that one
    stationary hit-and-run update in direction \(\theta\) changes the
    Boolean value.  Thus a lower bound for its spherical average is a
    Boolean conductance estimate for the hit-and-run kernel.

2.  The proposed bound

    \[
                       \sqrt n\,\mathbb E_\theta I_\theta\ge c             \tag{0.2}
    \]

    is false even when \(E\) is an exact attained balanced Cheeger
    minimizer.  For the isotropic cube and its coordinate half-cube,

    \[
       {1\over8n}\le\mathbb E_\theta I_\theta
          \le {1\over n}.                                    \tag{0.3}
    \]

    The coordinate half-cube is a half-volume isoperimetric minimizer, and
    concavity and symmetry of the convex-body isoperimetric profile make it
    a Cheeger minimizer.  Short random cube chords cause the extra
    \(n^{-1/2}\) loss.  A regular-simplex barycentric cap independently
    obeys \(1/[4(n+1)]\le\mathbb E I_\theta\le1/(n+1)\).

3.  Crofton plus the one-dimensional log-concave Cheeger inequality gives
    an exact scale-weighted quantity.  If \(\sigma_\theta(y)^2\) is the
    conditional variance along the line, then

    \[
      c_nP_\mu(E)\ge {1\over2\sqrt3}
       \mathbb E_{\theta,y}{\min(p_\theta,1-p_\theta)\over
                                  \sigma_\theta(y)},
      \qquad
      c_n=\mathbb E|\Theta_1|\asymp n^{-1/2}.                 \tag{0.4}
    \]

    Isotropy alone converts the right side only to

    \[
       P_\mu(E)\ge\sqrt{n/6}\,
                   \{\mathbb E_\theta I_\theta\}^{3/2}.      \tag{0.5}
    \]

    Hence even if (0.2) were true, this chain would yield only
    \(P_\mu(E)\ge c n^{-1/4}\).  The missing assertion is the weighted
    mixing bound

    \[
       \mathbb E_{\theta,y}{\min(p_\theta,1-p_\theta)\over
                                  \sigma_\theta(y)}
                              \ge {c\over\sqrt n},             \tag{0.6}
    \]

    not an unweighted hit-and-run conductance estimate.

4.  Beckmann transport does not bypass this loss.  If
    \(\mu_+=2\mu|_E\) and \(\mu_-=2\mu|_{E^c}\), then

    \[
       W_1(\mu_+,\mu_-)
       =\inf\left\{\int|b|\,d\mu:
            -\operatorname {div}(\mu b)=2q\mu\right\}.       \tag{0.7}
    \]

    The balanced long-tube estimate gives
    \(W_1(\mu_+,\mu_-)\ge c/\psi\), whereas isotropy gives only the
    upper bound \(2\sqrt n\).  A universal upper bound for (0.7), uniformly
    over balanced \(E\), is equivalent up to constants to target (T3).
    It may not be assumed.

For an exact Cheeger minimizer, extremality adds global BV
quasiminimality, CMC/Jacobi stability on regular strata, and the Young
contact law.  The cube calculation proves that these facts do not imply
(0.2).  The simplex cap violates the Young law at its support edges, the
product-exponential maximum is improved by a simultaneous ridge bevel, and
tribes-like sets pay large perimeter.  Turning the remaining weighted
statement into a theorem requires a new
**fiber-majority completion-or-saving inverse**: small weighted mixing must
either force a one-dimensional marginal witness or produce a finite
volume-preserving competitor.  This is the Boolean version of the global
completion gap in the normal-tube route.

## 1. Conditional mixing is hit-and-run Boolean conductance

Fix \(\theta\in S^{n-1}\).  Write points as \(x=y+t\theta\) and disintegrate

\[
 d\mu(x)=d\bar\mu_\theta(y)\,d\nu_{\theta,y}(t).              \tag{1.1}
\]

Every \(\nu_{\theta,y}\) is a one-dimensional log-concave probability on
an interval.  Given \(X=y+t\theta\), a stationary hit-and-run update in
direction \(\theta\) keeps \(y\) and independently resamples
\(t'\sim\nu_{\theta,y}\).  Conditional on \(y\), the probability that the
two Boolean values differ is

\[
                              2p_\theta(y)(1-p_\theta(y)).     \tag{1.2}
\]

If \(a=\min(p,1-p)\), then

\[
                            a\le2p(1-p)\le2a.                 \tag{1.3}
\]

Therefore, if \(K_\theta\) is the fixed-direction hit-and-run kernel,

\[
 I_\theta\le
 \mathbb P\{q(X')\ne q(X)\mid\theta\}\le2I_\theta,          \tag{1.4}
\]

and

\[
 \langle q,(I-K_\theta)q\rangle_{L^2(\mu)}
 =4\mathbb E_y p_\theta(y)(1-p_\theta(y)).                   \tag{1.5}
\]

Thus (0.2) is precisely a Boolean expansion estimate of order
\(n^{-1/2}\) for the direction-averaged hit-and-run kernel.

There is also a deterministic fiber-majority approximation.  Define

\[
 C_\theta=\{y+t\theta:p_\theta(y)\ge1/2\}.                   \tag{1.6}
\]

Then

\[
                             \mu(E\mathbin\triangle C_\theta)=I_\theta.
                                                                    \tag{1.7}
\]

Small \(I_\theta\) therefore says that \(E\) is close in mass to a set
which is constant on \(\theta\)-fibers.  It gives no bound for
\(P(C_\theta)\): the majority function of \(y\) can oscillate arbitrarily
fast.  This is the first place where Boolean mixing alone loses the
physical boundary.

## 2. Crofton and the conditional scale

Let \(E\) have finite weighted perimeter and reduced-boundary normal \(N\).
Define the directional boundary integral

\[
                    J_\theta(E)=\int_{\partial^*E}
                              |N\cdot\theta|\,d\sigma_\mu.    \tag{2.1}
\]

The one-dimensional BV slicing theorem gives the exact identity

\[
 J_\theta(E)=\mathbb E_y
                   P_{\nu_{\theta,y}}(E\cap(y+\mathbb R\theta)).
                                                                    \tag{2.2}
\]

For a one-dimensional log-concave probability \(\nu\) of variance
\(\sigma^2\),

\[
 P_\nu(A)\ge{1\over2\sqrt3\,\sigma}
                   \min(\nu(A),1-\nu(A)).                    \tag{2.3}
\]

When \(\sigma=0\), the conditional Boolean value is constant and both
sides are interpreted as zero.  Consequently

\[
 J_\theta(E)\ge {1\over2\sqrt3}
       B_\theta,\qquad
 B_\theta:=\mathbb E_y{a_\theta(y)\over\sigma_\theta(y)},
 \quad a_\theta=\min(p_\theta,1-p_\theta).                   \tag{2.4}
\]

Averaging \(|N\cdot\theta|\) over the sphere gives

\[
 \mathbb E_\theta J_\theta(E)=c_nP_\mu(E),\qquad
 c_n={\Gamma(n/2)\over\sqrt\pi\,\Gamma((n+1)/2)},            \tag{2.5}
\]

and \(c_n\le n^{-1/2}\).

The only scale control supplied directly by isotropy is

\[
 \mathbb E_y\sigma_\theta(y)^2
 =\mathbb E\operatorname {Var}(\langle X,\theta\rangle
               \mid P_{\theta^\perp}X)\le1.                 \tag{2.6}
\]

Cauchy--Schwarz, \(0\le a_\theta\le1/2\), and (2.6) give

\[
\begin{split}
 I_\theta^2
 &\le B_\theta\,\mathbb E_y[a_\theta\sigma_\theta],\\
 \mathbb E_y[a_\theta\sigma_\theta]
 &\le\sqrt{\mathbb E a_\theta\,
                 \mathbb E(a_\theta\sigma_\theta^2)}
 \le\sqrt{I_\theta/2}.
\end{split}                                                  \tag{2.7}
\]

Thus

\[
                              B_\theta\ge\sqrt2\,I_\theta^{3/2}.  \tag{2.8}
\]

Combine (2.4), (2.5), Jensen, and (2.8):

\[
 c_nP_\mu(E)\ge{1\over\sqrt6}
                  \mathbb E_\theta I_\theta^{3/2}
 \ge{1\over\sqrt6}\{\mathbb E_\theta I_\theta\}^{3/2}.    \tag{2.9}
\]

Using \(c_n\le n^{-1/2}\) proves (0.5).  The exponent \(3/2\) records a
real correlation obstruction: all the mixed fibers may be precisely those
with large conditional variance.  Conversely, very short fibers can make
\(I_\theta\) much smaller than the scale-normalized quantity \(B_\theta\),
as the simplex calculation shows next.

## 3. Exact model calculations

### 3.1 Gaussian halfspace: the \(n^{-1/2}\) scale

Let \(\mu=\gamma_n\) and \(E=\{x_1\ge0\}\).  Put
\(a=|\theta_1|\).  Conditional on \(Y=P_{\theta^\perp}X\), the remaining
coordinate along \(\theta\) is an independent standard Gaussian.  A
two-dimensional Gaussian wedge calculation gives

\[
 \boxed{
 I_\theta={1\over\pi}\arcsin|\theta_1|.}                    \tag{3.1}
\]

Indeed, after standardizing the component of \(e_1\) in
\(\theta^\perp\), the minority probability is
\(\Phi(-\sqrt{1-a^2}|Z|/a)\).  Its expectation is the angular mass of a
wedge of angle \(2\arcsin a\), namely (3.1).  Therefore

\[
 {1\over\pi}\mathbb E|\Theta_1|
 \le\mathbb E_\theta I_\theta
 \le {1\over2}\mathbb E|\Theta_1|,                          \tag{3.2}
\]

so \(\mathbb E I_\theta\asymp n^{-1/2}\).  Here every conditional
variance equals one, and unweighted and weighted mixing have the same
scale.  The Gaussian halfspace is the model for which (0.2) is correctly
normalized.

### 3.2 The half-cube: an exact Cheeger-minimizer counterexample

For \(n\ge2\), let \(\mu\) be uniform on the isotropic cube

\[
                              K=[-\sqrt3,\sqrt3]^n
\]

and let \(E=\{x_1\ge0\}\).  The coordinate half-cube is an isoperimetric
minimizer at volume \(1/2\).  This follows from the standard reflection to
the flat torus and the Hadwiger--Barthe--Maurey contraction argument.
Symmetry and concavity of the isoperimetric profile of a convex body make
\(I(v)/v\) nonincreasing on \((0,1/2]\).  Hence the same half-cube is an
attained Cheeger minimizer, with

\[
                 P_\mu(E)={1\over2\sqrt3},\qquad
                 \psi_\mu={1\over\sqrt3}.                     \tag{3.3a}
\]

We compute its random-fiber mixing without using minimality.  Write
\(a=\sqrt3\).  Fix \(\theta\in S^{n-1}\).  At a boundary point
\(x=(0,x_2,\ldots,x_n)\), let \(t_+\) and \(t_-\) be the distances along
the line to the two cube endpoints.  Since the conditional law is uniform
on the chord,

\[
 I_\theta={|\theta_1|\over2a}\,
       \mathbb E_{x_2,\ldots,x_n}\min(t_+,t_-).                \tag{3.3b}
\]

Put \(D_j=a-|x_j|\) for \(j\ge2\).  These are independent uniform random
variables on \([0,a]\), and

\[
 \min(t_+,t_-)=\min\left\{{a\over|\theta_1|},
             \min_{2\le j\le n}{D_j\over|\theta_j|}\right\}.
                                                                    \tag{3.3c}
\]

Let \(A=\sum_{j=2}^n|\theta_j|\).  If
\(Z=\min_{j\ge2}D_j/|\theta_j|\), then

\[
 \mathbb P\{Z>t\}=\prod_{j=2}^n
          \left(1-{t|\theta_j|\over a}\right)_+.             \tag{3.3d}
\]

The upper bound \(\prod(1-u_j)\le e^{-\sum u_j}\) gives

\[
                 \mathbb E\min(a/|\theta_1|,Z)\le {a\over A}.
                                                                    \tag{3.3e}
\]

For the lower bound, put

\[
                 t_0=\min\left\{{a\over2A},{a\over|\theta_1|}\right\}.
\]

For \(0\le t\le t_0\),
\(\prod(1-u_j)\ge1-\sum u_j\ge1/2\), and therefore

\[
 \mathbb E\min(a/|\theta_1|,Z)\ge {t_0\over2}.               \tag{3.3f}
\]

Equations (3.3b), (3.3e), and (3.3f) imply

\[
 {1\over8}\min\left\{{|\theta_1|\over A},2\right\}
 \le I_\theta
 \le {1\over2}\min\left\{{|\theta_1|\over A},1\right\}.  \tag{3.3g}
\]

Let

\[
                 w_1={|\theta_1|\over\sum_{j=1}^n|\theta_j|}.
\]

Exchangeability gives \(\mathbb Ew_1=1/n\).  Since
\(|\theta_1|/A=w_1/(1-w_1)\),

\[
 \min\left\{{|\theta_1|\over A},2\right\}\ge w_1,
 \qquad
 \min\left\{{|\theta_1|\over A},1\right\}\le2w_1.         \tag{3.3h}
\]

Averaging (3.3g) proves (0.3).  Thus (0.2) fails in the literal extremal
generality in which it was proposed.

### 3.3 A regular-simplex cap: an independent \(n^{-1}\) model

Work first in the regular simplex

\[
 \Delta_n=\{x=(x_0,\ldots,x_n):x_i\ge0,\ \sum_{i=0}^nx_i=1\}
\]

with its uniform probability on the affine hyperplane
\(H=\{\sum x_i=1\}\).  Its covariance on the tangent space
\(H_0=\{\sum z_i=0\}\) is a scalar matrix.  Passing to isotropic position
is therefore a scalar dilation, which leaves \(I_\theta\) and Haar-random
directions unchanged.

Let

\[
 E_s=\{x_0\ge s\},\qquad (1-s)^n={1\over2}.                  \tag{3.3}
\]

Fix a unit \(\theta\in H_0\).  On a line which crosses the interface
\(\Sigma_s=\{x_0=s\}\), let \(t_+\) and \(t_-\) be the distances from the
crossing point to the two endpoints of the simplex chord.  Since the
conditional law is uniform on that chord, its contribution to
\(I_\theta\), after multiplication by the fiber marginal, is
\(\min(t_+,t_-)\).  Projection of \(\Sigma_s\) to
\(\theta^\perp\cap H\) therefore gives

\[
 I_\theta=f_{X_0}(s)|\theta_0|\,
   \mathbb E_{x\in\Sigma_s}
       \min_{0\le i\le n}{x_i\over|\theta_i|},               \tag{3.4}
\]

where zero denominators are ignored and

\[
                         f_{X_0}(s)=n(1-s)^{n-1}.              \tag{3.5}
\]

To verify (3.4), note that the unit normal to \(\Sigma_s\) in \(H\) is
the normalization of the tangent gradient of \(x_0\).  The coarea factor
and the projection Jacobian multiply to
\(f_{X_0}(s)|\theta_0|\).  Also

\[
 \min(t_+,t_-)=\min_i{x_i\over|\theta_i|}.                   \tag{3.6}
\]

Put

\[
             A=\sum_{i=1}^n|\theta_i|,\qquad
             r={|\theta_0|\over A}.                           \tag{3.7}
\]

Conditional on \(x_0=s\), the remaining coordinates are uniform on the
simplex of total mass \(1-s\).  Hence, if
\(Z=\min_{1\le i\le n}x_i/|\theta_i|\),

\[
 \mathbb P\{Z>t\}
 =\left(1-{tA\over1-s}\right)^{n-1},\qquad
 0\le t\le{1-s\over A},                                    \tag{3.8}
\]

and

\[
                         \mathbb EZ={1-s\over nA}.             \tag{3.9}
\]

The additional cap \(s/|\theta_0|\) changes (3.9) by the factor

\[
 F_\theta=1-\left(1-\min\left\{1,
             {s\over(1-s)r}\right\}\right)^n.               \tag{3.10}
\]

Because \(\sum_i\theta_i=0\), one has \(0<r\le1\).  Moreover

\[
 {s\over1-s}=2^{1/n}-1\ge{\log2\over n},                    \tag{3.11}
\]

so (3.10) and \((1-u)^n\le e^{-nu}\) give

\[
                              {1\over2}\le F_\theta\le1.     \tag{3.12}
\]

Substituting (3.5), (3.9), and \((1-s)^n=1/2\) in (3.4) yields

\[
                              {r\over4}\le I_\theta\le{r\over2}.  \tag{3.13}
\]

Finally put

\[
 w_i={|\theta_i|\over\sum_{j=0}^n|\theta_j|}.
\]

Exchangeability gives \(\mathbb Ew_0=1/(n+1)\).  The zero-sum constraint
gives \(w_0\le1/2\), and

\[
                  w_0\le r={w_0\over1-w_0}\le2w_0.          \tag{3.14}
\]

Averaging (3.13) proves (0.3).

This is not a vanishing-perimeter construction.  In isotropic position the
same cap has

\[
 {P(E_s)\over1/2}
 ={n\sqrt n\over(n+1)\sqrt{n+2}},2^{1/n}\longrightarrow1.  \tag{3.15}
\]

It refutes (0.2) because the latter was proposed as a Boolean/fiber fact.
The flat cap is not asserted to be the exact relative isoperimetric region
of the simplex: it fails the orthogonal Young contact condition at its
support edges.  Thus exact Cheeger extremality would have to enter through
a new global contact/competitor theorem, not through Boolean balance or
isotropy.

### 3.4 Product-exponential maximum

For independent unit exponentials and the balanced set
\(E=\{\max_iY_i\ge L\}\),

\[
 P(E)={n\over2}(2^{1/n}-1)\in[\log2/2,1/2],\qquad
 Q={I_n\over n}.                                             \tag{3.16}
\]

Crofton therefore gives the exact aggregate scale

\[
                       \mathbb E_\theta J_\theta(E)=c_nP(E)\asymp n^{-1/2}.
                                                                    \tag{3.17}
\]

The conditional fibers are truncated by the orthant walls and by competing
maximum facets; their conditional variances cannot be replaced by one as
in the Gaussian calculation.  The already audited normal-ray calculation
shows the same phenomenon geometrically: short killed loss is small while
the full-slice completion defect tends to the entire perimeter.  Thus
(3.16) is compatible with the weighted target (0.6), but it invalidates
any step which silently identifies \(I_\theta\) with \(J_\theta\).

This maximum set is not an exact minimizer.  A simultaneous fixed ridge
bevel lowers its half-mass perimeter by more than `.0061` for large \(n\).
The Boolean inverse needed for a true minimizer must detect that finite
ridge operation; the regular fiber conditionals do not.

### 3.5 Gaussian tribes

Let the signs of independent standard Gaussian coordinates be unbiased
bits and let \(E\) be the usual balanced tribes Boolean function: an OR of
\(m\) disjoint AND-blocks of size \(b\), with
\(m\asymp2^b\) and \(n=mb\).  Its total discrete influence is
\(\Theta(b)=\Theta(\log n)\).  Each pivotal coordinate contributes the
Gaussian density \((2\pi)^{-1/2}\) at zero, so

\[
                              P_{\gamma_n}(E)=\Theta(\log n). \tag{3.18}
\]

Every Gaussian conditional line has variance one.  Equations (2.2)--(2.4)
therefore imply

\[
             \mathbb E_\theta I_\theta
             \le2\sqrt3\,c_nP_{\gamma_n}(E)
             =O\left({\log n\over\sqrt n}\right).            \tag{3.19}
\]

Tribes is genuinely nonlinear and has more boundary/mixing than a
halfspace by the expected logarithmic factor.  It is excluded from
Cheeger extremality by its large perimeter, not by a general theorem saying
that every Boolean function has halfspace-like fiber mixing.

## 4. What exact Cheeger extremality adds

If \(E\) is an attained balanced Cheeger minimizer with
\(P(E)=\psi/2\), then for every finite-perimeter \(F\)

\[
             P(E)\le P(F)+\psi\,\mu(E\mathbin\triangle F),    \tag{4.1}
\]

and \(P(E)\le P(F)\) whenever \(\mu(F)=1/2\).  On smooth regular
strata this gives constant weighted mean curvature, volume-constrained
Jacobi stability, and the natural Young contact condition.

Apply (4.1) to the fiber-majority set (1.6):

\[
                  P(E)\le P(C_\theta)+\psi I_\theta.          \tag{4.2}
\]

Equation (4.2) is exact but has the wrong direction unless one can control
the transverse perimeter of \(C_\theta\).  Neither the hit-and-run flip
probability nor conditional variance controls that perimeter.  The simplex
calculation makes the issue quantitative: a random line changes the Boolean
value with probability \(\Theta(1/n)\), although the physical perimeter is
\(\Theta(1)\).

A sufficient new inverse would be the following.

> **Fiber-majority completion-or-saving theorem (missing).**  For every
> balanced Cheeger minimizer, either
> \[
>    \mathbb E_{\theta,y}{a_\theta(y)\over\sigma_\theta(y)}
>       \ge {c_0\over\sqrt n},                                \tag{4.3}
> \]
> or there is a direction \(u\), a median \(t\) of
> \(\langle X,u\rangle\), and a halfspace
> \(H=\{\langle x,u\rangle\ge t\}\) such that the finite replacement of
> \(E\) by \(H\), with exact mass correction, has strictly smaller
> perimeter.  In a packet formulation, large transverse perimeter of the
> majority cylinders must be converted with bounded reuse into a
> ridge/contact/medial saving.

For a minimizer the saving alternative is impossible; (4.3), (2.4), and
(2.5) give a universal lower bound for \(P(E)\).  Proving (4.3) directly is
already a dimension-free isoperimetric statement.  Calling nearly pure
fibers a product, or discarding their transverse majority boundary, does
not prove it.

## 5. Beckmann and Wasserstein-1

Put

\[
                         \mu_+=2\,1_E\mu,\qquad
                         \mu_-=2\,1_{E^c}\mu.                 \tag{5.1}
\]

Kantorovich duality and the Beckmann formulation give

\[
\begin{split}
 W_E:=W_1(\mu_+,\mu_-)
 &=2\sup_{\operatorname {Lip}(f)\le1}\int qf\,d\mu\\
 &=\inf\left\{\int|b|\,d\mu:
       -\operatorname {div}(\mu b)=2q\mu\right\},
\end{split}                                                  \tag{5.2}
\]

where the second line is interpreted in distributions and may equivalently
use vector-valued finite measures.

Let \(d_E\) be signed distance from the interface, positive on \(E\).  It
is one-Lipschitz and

\[
                              W_E\ge2\int|d_E|\,d\mu.          \tag{5.3}
\]

For a balanced exact minimizer, apply the short two-sided tube estimate with
any fixed \(0<\gamma<1/4\).  The distance from the interface exceeds
\(T\ge\log(1+\gamma)/\psi\) on mass at least
\(1-\gamma\).  Therefore

\[
                    W_E\ge {2(1-\gamma)\log(1+\gamma)\over\psi}.
                                                                    \tag{5.4}
\]

On the other hand, the independent coupling and isotropy give only

\[
 W_E\le\mathbb E|X_+-X_-|
 \le\mathbb E|X_+|+\mathbb E|X_-|
 =2\mathbb E|X|\le2\sqrt n.                                 \tag{5.5}
\]

This recovers only \(\psi\ge c/\sqrt n\).

The desired dimension-free upper bound on \(W_E\) is not an independent
transport lemma.  If target (T3) holds with constant \(C\), then for every
one-Lipschitz \(f\), choosing a median \(m_f\) gives

\[
 \left|2\int qf\,d\mu\right|
 =\left|2\int q(f-m_f)\,d\mu\right|
 \le2\int|f-m_f|\,d\mu\le2C,                                \tag{5.6}
\]

so \(W_E\le2C\).  Conversely, suppose \(W_E\le C\) for every balanced
Borel \(E\).  For a one-Lipschitz \(f\) with a non-atomic median, take
\(E=\{f\ge m_f\}\).  Then

\[
 2\int|f-m_f|\,d\mu
 =\int f\,d(\mu_+-\mu_-)\le W_E\le C.                       \tag{5.7}
\]

Level splitting handles a median atom.  Thus uniform balanced-set
Beckmann control is equivalent, up to a factor two, to first-moment
concentration.  It cannot be used as an input to prove KLS.

## 6. Audit conclusion

The unweighted Boolean estimate (0.2) is disproved, even for an exact
balanced Cheeger minimizer, by the half-cube calculation (0.3).  The simplex
gives an independent sharp short-fiber model.  Even if (0.2) were postulated
for some narrower class, the only unconditional
fiber-scale estimate loses the exponent in (0.5).  Beckmann transport
repackages target (T3) rather than proving it.

The viable Boolean target is the scale-weighted inverse (4.3), together
with a proof that failure of it creates a genuine half-mass perimeter
competitor.  The required competitor must retain support contacts, fiber
completion, and transverse majority perimeter.  These are exactly the
quantities missed by unweighted hit-and-run mixing.
