# Weighted random fibers for a balanced Cheeger minimizer

## Outcome

Let \(\mu\) be a full-dimensional isotropic log-concave probability on
\(\mathbb R^n\), with density \(\rho=e^{-V}\), and let \(E\) be a
finite-perimeter set with \(\mu(E)=1/2\).  For
\(\theta\in S^{n-1}\), write \(x=y+t\theta\), \(y\in\theta^\perp\), and
put

\[
\begin{split}
 w_\theta(y)&=\int_{\mathbb R}\rho(y+t\theta)\,dt,\\
 m_{\theta,y}^{+}
   &=\int_{\mathbb R}1_E(y+t\theta)\rho(y+t\theta)\,dt,\\
 m_{\theta,y}^{-}&=w_\theta(y)-m_{\theta,y}^{+}.
\end{split}                                                    \tag{0.1}
\]

If \(\nu_{\theta,y}\) is the conditional probability on the line and
\(\sigma_\theta(y)\) its standard deviation, define

\[
 B_\theta(E)=\int_{\theta^\perp}
   {\min(m_{\theta,y}^{+},m_{\theta,y}^{-})
       \over \sigma_\theta(y)}\,dy.                            \tag{0.2}
\]

The normalization in (0.2) is essential: the two masses are
*unnormalized*, whereas \(\sigma_\theta(y)\) belongs to the normalized
conditional probability.  Equivalently,

\[
 B_\theta(E)=\int_{\theta^\perp}
       w_\theta(y){\min(p_\theta(y),1-p_\theta(y))
       \over\sigma_\theta(y)}\,dy,                             \tag{0.3}
\]

where \(p_\theta(y)=m_{\theta,y}^{+}/w_\theta(y)\).

The exact slicing and Crofton calculation gives

\[
 J_\theta(E):=\int_{\partial^*E}|N_E\!\cdot\theta|\,d\sigma_\mu
 \ge {1\over2\sqrt3}\,B_\theta(E),                            \tag{0.4}
\]

and

\[
 \mathbb E_{\theta\sim{\rm unif}(S^{n-1})}J_\theta(E)
   =c_nP_\mu(E),\qquad
 c_n={\Gamma(n/2)\over\sqrt\pi\,\Gamma((n+1)/2)},\qquad
 {1\over\sqrt{2n}}\le c_n\le {1\over\sqrt n}.                \tag{0.5}
\]

Consequently, the proposed estimate

\[
             \mathbb E_\theta B_\theta(E)\ge{\beta\over\sqrt n}
                                                                  \tag{WFC}
\]

for every balanced exact Cheeger minimizer would imply

\[
 P_\mu(E)\ge{\beta\over2\sqrt3},\qquad
 \psi_\mu=2P_\mu(E)\ge{\beta\over\sqrt3}.                     \tag{0.6}
\]

Thus (WFC) is a genuinely KLS-closing statement, not an auxiliary
normalization estimate.

The principal new conclusions of this report are:

1.  (0.4) is proved below with the explicit constant and all marginal
    factors visible.
2.  Exact Cheeger minimality permits fiberwise mass-preserving
    replacement, but only controls the *total* perimeter.  The
    one-dimensional replacement controls \(J_\theta\), while its
    uncontrolled transverse graph term is exactly the missing quantity.
3.  Majority completion creates a mass defect at most the unweighted
    fiber error.  The projected max-flow calibration gives exactly
    \(P_{\bar\mu_\theta}(\{p_\theta>1/2\})
      \ge P_\mu(E)-\psi_\mu I_\theta\), so naive dimension descent merely
    reproduces quasiminimality.  A strict repair needs new transverse
    information.
4.  The target survives the two available exact-minimizer stress tests.
    For a Gaussian halfspace,

    \[
       B_\theta={1\over\pi}\arcsin|\theta_1|;
                                                                  \tag{0.7}
    \]

    for the coordinate half-cube in the isotropic cube,

    \[
       {\sqrt{12}\over192}J_\theta
          \le B_\theta\le\sqrt3\,J_\theta,\qquad
       {c_n\over192}\le\mathbb E_\theta B_\theta\le{c_n\over2}.
                                                                  \tag{0.8}
    \]

    Hence neither is a counterexample.
5.  The regular-simplex median cap obeys the same comparison as (0.8),
    with \(192\) unchanged.  The radial-exponential median sphere also
    has \(\mathbb E_\theta B_\theta\asymp n^{-1/2}\).  These are useful
    stress models but are not invoked as exact global minimizers.  The
    product-exponential maximum is fiber-efficient in every coordinate
    direction but fails exact minimality by a ridge bevel.

No proof of (WFC), and no exact-minimizer counterexample to it, is obtained.
The remaining load-bearing statement is isolated in Section 7 as a
finite, quantitative completion-or-saving theorem.  Its saving alternative
must control the transverse perimeter without assuming KLS for a marginal.

## 1. Disintegration and normalization

### 1.1 Conditional probabilities

Fix \(\theta\in S^{n-1}\).  Fubini gives

\[
 \int_{\mathbb R^n}g(x)\rho(x)\,dx
 =\int_{\theta^\perp}\int_{\mathbb R}
       g(y+t\theta)\rho(y+t\theta)\,dt\,dy.                   \tag{1.1}
\]

The marginal density of \(Y=P_{\theta^\perp}X\) is \(w_\theta(y)\).
On \(\{w_\theta>0\}\),

\[
 d\nu_{\theta,y}(t)
  ={\rho(y+t\theta)\over w_\theta(y)}\,dt.                   \tag{1.2}
\]

Prékopa's theorem says that \(w_\theta\) is log-concave on
\(\theta^\perp\), and restriction of \(V\) to a line says that every
\(\nu_{\theta,y}\) is a one-dimensional log-concave probability on an
interval.  Put

\[
 \bar t_\theta(y)=\int t\,d\nu_{\theta,y}(t),\qquad
 \sigma_\theta(y)^2
   =\int(t-\bar t_\theta(y))^2\,d\nu_{\theta,y}(t).           \tag{1.3}
\]

Isotropy and conditional variance decomposition give

\[
 \int_{\theta^\perp}w_\theta(y)\sigma_\theta(y)^2\,dy
 =\mathbb E\operatorname {Var}(\langle X,\theta\rangle\mid Y)
 \le1.                                                       \tag{1.4}
\]

For a full-dimensional density, \(\sigma_\theta(y)>0\) for
\(w_\theta(y)\,dy\)-almost every \(y\) in the relative interior of the
marginal support.  On any exceptional zero-variance fiber the conditional
Boolean value is constant; its contribution to (0.2) is defined to be
zero.

### 1.2 Why the marginal factor cannot be omitted

Let \(q_y(t)=\rho(y+t\theta)\), \(w=\int q_y\), and let
\(\widetilde q_y=cq_y\) for a positive constant \(c\).  The conditional
probability and its standard deviation do not change, while both
unnormalized phase masses and the weighted one-dimensional perimeter are
multiplied by \(c\).  Thus the scale-covariant fiber quantity is

\[
 {\min(\int_{E_y}q_y,\int_{E_y^c}q_y)\over\sigma_y}
   ={w_y\min(\nu_y(E_y),1-\nu_y(E_y))\over\sigma_y}.          \tag{1.5}
\]

Using \(\min(p_y,1-p_y)/\sigma_y\) against Lebesgue \(dy\), without
\(w_y\), is not invariant under this elementary rescaling and cannot be
compared to \(P_\mu(E)\).

For later use define the unweighted conditional minority mass

\[
 I_\theta(E)=\int_{\theta^\perp}
          w_\theta(y)a_\theta(y)\,dy,\qquad
 a_\theta(y)=\min(p_\theta(y),1-p_\theta(y)).                 \tag{1.6}
\]

From (1.4), Cauchy--Schwarz, and \(a_\theta\le1/2\),

\[
\begin{split}
 I_\theta^2
 &\le B_\theta
       \int w_\theta a_\theta\sigma_\theta,\\
 \int w_\theta a_\theta\sigma_\theta
 &\le\left(I_\theta
       \int w_\theta a_\theta\sigma_\theta^2\right)^{1/2}
 \le\sqrt{I_\theta/2}.
\end{split}                                                   \tag{1.7}
\]

Therefore

\[
                     B_\theta\ge\sqrt2\,I_\theta^{3/2},
 \qquad
 I_\theta\le2^{-1/3}B_\theta^{2/3}.                          \tag{1.8}
\]

The exponent \(3/2\) cannot be discarded by covariance alone: in the
half-cube, \(\mathbb E_\theta I_\theta\asymp n^{-1}\) but
\(\mathbb E_\theta B_\theta\asymp n^{-1/2}\).

## 2. The one-dimensional constant

### 2.1 Median density

**Lemma 2.1.**  Let \(\nu\) be a nondegenerate log-concave probability on
\(\mathbb R\), with density \(f\), median \(m\), and standard deviation
\(\sigma\).  Then

\[
                  h_\nu:=\inf_A
     {P_\nu(A)\over\min(\nu(A),1-\nu(A))}
     =2f(m)\ge{1\over2\sqrt3\,\sigma}.                       \tag{2.1}
\]

Here the infimum is over Borel sets of nontrivial mass and \(P_\nu\) is
the exterior Minkowski content, equivalently the weighted BV perimeter.

**Proof.**  The one-dimensional log-concave isoperimetric theorem says
that, at prescribed mass, a half-line minimizes weighted perimeter.  It
can be seen directly by ordering the boundary points of a finite union of
intervals and using log-concavity; approximation then gives the Borel
statement.

The distribution function \(F\) and the survival function \(1-F\) are
log-concave.  Hence \(f/F=(\log F)'\) is nonincreasing and
\(f/(1-F)=-(\log(1-F))'\) is nondecreasing.  For \(x\le m\) and
\(x\ge m\), respectively,

\[
 {f(x)\over F(x)}\ge2f(m),\qquad
 {f(x)\over1-F(x)}\ge2f(m).                                 \tag{2.2}
\]

Thus every half-line has perimeter-to-minority ratio at least \(2f(m)\),
and a median half-line attains equality.

A standard one-dimensional consequence of unimodality and log-concavity
is

\[
                         \|f\|_\infty\le2f(m).               \tag{2.3}
\]

Here is a direct proof.  Normalize a mode to \(x_0=0\) and
\(\|f\|_\infty=f(0)=1\), and suppose \(0<m\); reflection handles the
other case.  Put \(a=f(m)\) and \(s=-\log(a)/m\).  Concavity of
\(\log f\) gives

\[
 f(x)\ge e^{-sx}\quad(0\le x\le m),\qquad
 f(m+t)\le ae^{-st}\quad(t\ge0).                             \tag{2.3a}
\]

The second inequality follows because every right derivative of
\(\log f\) after \(m\) is at most the chord slope
\((\log f(m)-\log f(0))/m=-s\).  Since the two median half-lines have
equal mass,

\[
 {1-a\over s}\le\int_0^mf(x)\,dx
 \le\int_{-\infty}^mf(x)\,dx
 =\int_m^\infty f(x)\,dx
 \le{a\over s}.                                              \tag{2.3b}
\]

Thus \(a\ge1/2\), proving (2.3).  Plateau modes and \(m=0\) follow by a
limit and are easier.

Finally, if a probability density is bounded by \(M\), then for every
center \(c\),

\[
 \mathbb P\{|X-c|>r\}\ge(1-2Mr)_+.
\]

Integrating the tail of \((X-c)^2\) gives

\[
 \mathbb E(X-c)^2
 \ge\int_0^{1/(2M)}2r(1-2Mr)\,dr={1\over12M^2}.              \tag{2.3c}
\]

Taking \(c=\mathbb EX\) and \(M=\|f\|_\infty\), we obtain

\[
 \sigma^2\ge{1\over12\|f\|_\infty^2}
      \ge{1\over48f(m)^2}.                                   \tag{2.4}
\]

Equations (2.1)--(2.4) prove the claim. \(\square\)

The same statement for the unnormalized line density \(q_y=w_yf_y\) is

\[
 P_{q_y}(A)\ge {1\over2\sqrt3\,\sigma_y}
       \min\left(\int_Aq_y,\int_{A^c}q_y\right).              \tag{2.5}
\]

### 2.2 Slicing

Let \(E\) have finite \(\mu\)-weighted perimeter and let \(N_E\) be its
measure-theoretic exterior unit normal.  Define

\[
 J_\theta(E)=\int_{\partial^*E}|N_E\!\cdot\theta|
                       \rho\,d\mathcal H^{n-1}.              \tag{2.6}
\]

The BV slicing theorem applied to the scalar measure
\(\theta\cdot D1_E\), followed by Fubini with the locally integrable
weight \(\rho\), gives

\[
 J_\theta(E)=\int_{\theta^\perp}
      P_{q_{\theta,y}}(E_y)\,dy.                              \tag{2.7}
\]

For a general log-concave density, \(\rho\) is continuous on the relative
interior of its convex support.  Apply (2.7) first on compact subsets of
that interior and then increase to the support; the boundary contribution
is the usual relative perimeter.  No smoothness or strict convexity of
\(V\) is used.

Combining (2.5) and (2.7) proves (0.4).

### 2.3 Crofton normalization

Fubini and rotational invariance give

\[
\begin{split}
 \mathbb E_\theta J_\theta(E)
 &=\int_{\partial^*E}\mathbb E_\theta|N_E(x)\!\cdot\theta|
                      \,d\sigma_\mu(x)\\
 &=c_nP_\mu(E).
\end{split}                                                   \tag{2.8}
\]

If \(\Theta\) is uniform on \(S^{n-1}\), beta integration gives

\[
 c_n=\mathbb E|\Theta_1|
   ={\Gamma(n/2)\over\sqrt\pi\,\Gamma((n+1)/2)}.             \tag{2.9}
\]

The upper bound \(c_n\le n^{-1/2}\) follows immediately from
\(\mathbb E\Theta_1^2=1/n\).  Wendel's gamma-ratio inequality gives
\(c_n\ge\sqrt{2/(\pi n)}\ge(2n)^{-1/2}\).  This completes the exact
proof of (0.4)--(0.5).

## 3. What exact Cheeger minimality supplies

Assume now that \(E\) is an attained balanced Cheeger minimizer:

\[
          \mu(E)=\frac12,\qquad P_\mu(E)=\frac{\psi_\mu}{2}. \tag{3.1}
\]

For every finite-perimeter \(F\),

\[
                 P_\mu(E)\le
        P_\mu(F)+\psi_\mu\,\mu(E\mathbin\triangle F).         \tag{3.2}
\]

Indeed,

\[
\begin{split}
 P_\mu(F)
 &\ge\psi_\mu\min(\mu(F),1-\mu(F))\\
 &\ge\psi_\mu\left({1\over2}
          -|\mu(F)-1/2|\right)\\
 &\ge P_\mu(E)-\psi_\mu\mu(E\mathbin\triangle F).
\end{split}                                                   \tag{3.3}
\]

In particular, \(P_\mu(E)\le P_\mu(F)\) whenever \(\mu(F)=1/2\).
This is the exact global information available for replacement.

### 3.1 Fiberwise isoperimetric replacement

For fixed \(\theta\), choose on almost every fiber a half-line
\(H_y\subset\mathbb R\) with

\[
                    \nu_{\theta,y}(H_y)=p_\theta(y)          \tag{3.4}
\]

and with the smaller of the lower-tail and upper-tail perimeters.  Quantile
measurability makes

\[
 R_\theta E=\{y+t\theta:t\in H_y\}                           \tag{3.5}
\]

measurable.  It preserves every fiber mass, so
\(\mu(R_\theta E)=1/2\), and the one-dimensional isoperimetric theorem
gives

\[
                         J_\theta(R_\theta E)\le J_\theta(E).
                                                                  \tag{3.6}
\]

Equation (3.6) is not a total-perimeter comparison.  Even when
\(R_\theta E\) is an epigraph \(t>g(y)\), its weighted perimeter is

\[
 P_\mu(R_\theta E)
   =\int_{\theta^\perp}\rho(y+g(y)\theta)
       \sqrt{1+|\nabla g(y)|^2}\,dy,                          \tag{3.7}
\]

whereas its directional perimeter is only

\[
 J_\theta(R_\theta E)
   =\int_{\theta^\perp}\rho(y+g(y)\theta)\,dy.                \tag{3.8}
\]

Neither \(B_\theta\), \(I_\theta\), nor the conditional variances control
\(\nabla g\).  Exact minimality applied to (3.5) therefore reads

\[
 P_\mu(E)\le P_\mu(R_\theta E),                              \tag{3.9}
\]

but (3.6) supplies no upper bound for the right side.

This is not a technical defect of choosing quantiles.  A balanced
cylinder \(C=A+\mathbb R\theta\) has \(B_\theta(C)=0\) while

\[
                         P_\mu(C)=P_{\bar\mu_\theta}(A),      \tag{3.10}
\]

which may be positive or arbitrarily large.  Thus no directional estimate
of transverse perimeter by \(B_\theta\) is possible.

### 3.2 Majority completion and its mass defect

Define the majority cylinder

\[
 C_\theta=\{y+t\theta:p_\theta(y)\ge1/2\}.                   \tag{3.11}
\]

Then

\[
 \mu(E\mathbin\triangle C_\theta)=I_\theta,\qquad
 |\mu(C_\theta)-1/2|\le I_\theta.                            \tag{3.12}
\]

If \(C_\theta\) has finite perimeter, (3.2) gives

\[
 P_\mu(E)\le P_{\bar\mu_\theta}(\{p_\theta\ge1/2\})
                      +\psi_\mu I_\theta.                   \tag{3.13}
\]

This has the wrong direction unless the transverse majority perimeter is
controlled.  Moreover, making the cylinder exactly balanced requires
altering the base set by marginal mass at most \(I_\theta\).  The marginal
\(\bar\mu_\theta\) is log-concave and has covariance at most the identity,
but neither a Cheeger *lower* bound nor the standard halfspace *upper*
bound supplies a repair with the sign and size needed in (3.13).
Postulating a balanced marginal replacement whose perimeter is below
\(P_\mu(E)\) simply postulates the saving alternative in one lower
dimension.

### 3.3 The projected max-flow calibration is tautological

There is an exact dual explanation for the sign in (3.13).  Put

\[
                        q=1_E-1_{E^c}.                       \tag{3.13a}
\]

For every finite-perimeter \(A\),

\[
 \psi_\mu\left|\int_Aq\,d\mu\right|
 \le\psi_\mu\min(\mu(A),1-\mu(A))
 \le P_\mu(A).                                               \tag{3.13b}
\]

The continuous max-flow/min-cut duality therefore produces
\(z\in L^\infty(\mu;\mathbb R^n)\), \(|z|\le1\), such that

\[
                     -\operatorname {div}(\rho z)
                         =\psi_\mu q\rho                     \tag{3.13c}
\]

in distributions.  Here is a direct functional-analytic verification.
By layer cake, (3.13b), and weighted coarea,
\[
 \left|\psi_\mu\int qu\,d\mu\right|\le {\rm TV}_\mu(u)
\]
for every compactly supported weighted-BV function \(u\).  Hence the
functional \(\nabla u\mapsto\psi_\mu\int qu\,d\mu\) has norm at most one
on the subspace of gradients in \(L^1(\mu;\mathbb R^n)\).  Hahn--Banach
and \(L^1\)-\(L^\infty\) duality give \(z\) and (3.13c).  Since
\(\psi_\mu\int_Eq\,d\mu=\psi_\mu/2=P_\mu(E)\), the normal trace of \(z\)
also saturates capacity on \(\partial^*E\).

Project this flow onto \(\theta^\perp\):

\[
 \bar z_\theta(y)
 =\mathbb E[P_{\theta^\perp}z(X)\mid
                P_{\theta^\perp}X=y],\qquad
 r_\theta(y)=\mathbb E[q(X)\mid
                P_{\theta^\perp}X=y]=2p_\theta(y)-1.         \tag{3.13d}
\]

Then \(|\bar z_\theta|\le1\), and testing (3.13c) with functions of \(y\)
gives

\[
 -\operatorname {div}_{\theta^\perp}
       (w_\theta\bar z_\theta)
       =\psi_\mu w_\theta r_\theta.                          \tag{3.13e}
\]

Let \(A_\theta=\{r_\theta>0\}\).  If it has finite marginal perimeter,
test (3.13e) with \(2\,1_{A_\theta}-1\).  The capacity bound yields

\[
\begin{split}
 2P_{\bar\mu_\theta}(A_\theta)
 &\ge\psi_\mu\int|r_\theta|w_\theta\,dy\\
 &=\psi_\mu(1-2I_\theta),
\end{split}                                                   \tag{3.13f}
\]

or

\[
 P_{\bar\mu_\theta}(A_\theta)
 \ge P_\mu(E)-\psi_\mu I_\theta.                             \tag{3.13g}
\]

This is exactly (3.13), rearranged.  Thus projecting the extremal
calibration onto the majority base does not create a dimension-descent
gain: it certifies the same error budget as quasiminimality.  Any strict
saving must exploit information discarded by the projection, such as
transverse coherence across several directions or the support-contact
flux of Section 3.5.

### 3.4 A precise gluing formula

Let \(A\subset\theta^\perp\) be a bounded Lipschitz set and let \(F\)
coincide with \(E\) outside the cylinder
\(\mathcal C=A+\mathbb R\theta\).  The BV gluing inequality is

\[
\begin{split}
 P_\mu(F)
 \le{}&P_\mu(E;\mathcal C^c)+P_\mu(F;\mathcal C)\\
 &+\int_{\partial A}\int_{\mathbb R}
       |1_{F}(y+t\theta)-1_E(y+t\theta)|
       \rho(y+t\theta)\,dt\,d\mathcal H^{n-2}(y).
\end{split}                                                   \tag{3.14}
\]

The last term is the side-wall trace cost.  Randomly translating a slab
partition of mesh \(R\) in one base direction makes its expected side
cost \(2\mu(E\mathbin\triangle F)/R\).  This can make *external* gluing
cheap.  It does not control the internal transverse graph term in (3.7).

If instead one uses a single horizontal cut in each base cell, the
internal graph term vanishes, but only the aggregate cell mass is
preserved.  The location of that cut and the mismatch on cell walls are
not controlled by \(B_\theta\).  Passing to larger cells converges to a
global median halfspace, whose perimeter is universally bounded above,
not below.  It cannot contradict a hypothetical minimizer of very small
perimeter.

These two alternatives expose the exact replacement obstruction:
fiberwise mass preservation creates transverse oscillation, while
transverse smoothing creates a mass-repair problem.

### 3.5 Support-killed fibers and the contact tensor

There is a second bookkeeping issue when \(\mu\) is uniform on a smooth
convex body \(K\).  Let \(\Sigma\) be the regular interior interface of a
relative minimizer, \(\Gamma=\partial\Sigma\subset\partial K\), \(H\) its
constant mean curvature, and \(\nu\) its conormal on \(\Gamma\).  To avoid
confusion with the fiber functional, denote the support contact tensor by

\[
\begin{split}
 Q&={1\over P}\int_\Sigma N\otimes N\,dA,\\
 X&={1\over P}\int_\Sigma(x-c)\otimes N\,dA,\\
 \mathcal T_K&={1\over P}\int_\Gamma(x-c)\otimes\nu\,dS.
\end{split}                                                   \tag{3.15}
\]

The tensor Minkowski identity proved in
work/tensor_minkowski_completion.md is

\[
                           I-Q=HX+\mathcal T_K.              \tag{3.16}
\]

Since \(\operatorname {tr}Q=1\),

\[
 \sqrt{n-2}\le\|I-Q\|_F
 \le |H|\left({1\over P}\int_\Sigma|x-c|^2\,dA\right)^{1/2}
       +\|\mathcal T_K\|_F.                                  \tag{3.17}
\]

Thus, if the surface-position moment is at most \(C_0n\),
\(|H|\le2P\), and \(P\) is small, a tensor of order \(\sqrt n\) must be
carried by support contact.

This is exactly the branch that a fiber replacement must retain.  A line
whose minority segment is terminated rapidly by \(\partial K\) may have
small unweighted mixing.  In relative perimeter its support endpoint costs
nothing, but after majority completion its trace reappears in the
side-wall term of (3.14).  Discarding such endpoint-terminated fibers
discards the flux represented by \(\mathcal T_K\).  A valid completion
theorem must therefore give one of two outputs:

1. enough two-sided interior fiber mass to prove (WFC); or
2. a bounded-reuse family of support-contact cells whose finite
   replacement saves perimeter.

For the coordinate half-cube, \(H=0\), \(Q=e_1\otimes e_1\), and
\(\mathcal T_K=I-e_1\otimes e_1\).  The contact tensor carries all
transverse directions, while Section 5 shows that the weighted random
fibers are nevertheless efficient.  Support contact is therefore a
structured completion mechanism, not an error term.

## 4. Exact stress test I: a Gaussian halfspace

Let \(\mu=\gamma_n\) and \(E=\{x_1\ge0\}\).  This is an exact balanced
isoperimetric, hence Cheeger, minimizer.  Fix \(\theta\), put
\(a=|\theta_1|\), and reverse \(\theta\) if necessary.  With independent
standard Gaussians \(T,Z\),

\[
 X_1=aT+\sqrt{1-a^2}\,Z.                                    \tag{4.1}
\]

Conditioning on \(P_{\theta^\perp}X\) fixes \(Z\), while \(T\) remains
standard Gaussian.  Hence \(\sigma_\theta(y)=1\), and the conditional
minority probability is the Bayes error for predicting the sign of
\(X_1\) from \(Z\).  The sign-disagreement formula for correlated centered
Gaussians gives

\[
 B_\theta=I_\theta
 =\mathbb P\{X_1Z<0\}
 ={1\over\pi}\arccos\sqrt{1-a^2}
 ={1\over\pi}\arcsin a.                                     \tag{4.2}
\]

Therefore

\[
 {c_n\over\pi}
 \le\mathbb E_\theta B_\theta
 \le{c_n\over2}.                                             \tag{4.3}
\]

The lower bound uses \(\arcsin a\ge a\), and the upper bound uses
\(\arcsin a\le\pi a/2\).

## 5. Exact stress test II: the half-cube

Let \(a=\sqrt3\), let \(\mu\) be uniform on
\(K=[-a,a]^n\), and put \(E=\{x_1\ge0\}\).  The coordinate half-cube is a
half-volume isoperimetric minimizer (Hadwiger's cube theorem, equivalently
the Barthe--Maurey contraction proof).  Concavity of the relative
isoperimetric profile implies that it is also a Cheeger minimizer, with

\[
             P_\mu(E)={1\over2a},\qquad \psi_\mu={1\over a}. \tag{5.1}
\]

Fix \(\theta\) with \(\theta_1\ne0\).  A line meeting the interface has a
unique crossing point
\(x=(0,x_2,\ldots,x_n)\).  Let \(\ell_+(x)\) and \(\ell_-(x)\) be its two
chord lengths from \(x\) to \(\partial K\).  Since the conditional law is
uniform on a chord,

\[
 \sigma_\theta={\ell_++\ell_-\over\sqrt{12}},\qquad
 a_\theta={\min(\ell_+,\ell_-)\over\ell_++\ell_-}.            \tag{5.2}
\]

Projection of the interface onto \(\theta^\perp\) has Jacobian
\(|\theta_1|\).  Thus

\[
 B_\theta
 =\sqrt{12}\,J_\theta
   \mathbb E_{x_2,\ldots,x_n}
       {\min(\ell_+,\ell_-)\over\ell_++\ell_-},\qquad
 J_\theta={|\theta_1|\over2a}.                               \tag{5.3}
\]

We give a uniform lower bound for the expectation.  Put

\[
 L=|\theta_1|+\sum_{j=2}^n|\theta_j|,\qquad
 r_0={a\over4L}.                                             \tag{5.4}
\]

For \(j\ge2\), the distance from a uniform \(x_j\in[-a,a]\) to either
selected facet is uniform on \([0,2a]\).  A union bound gives

\[
 \mathbb P\{\ell_+\ge r_0,\ \ell_-\ge r_0\}
 \ge1-{r_0\over a}\sum_{j=2}^n|\theta_j|
 \ge{3\over4}.                                               \tag{5.5}
\]

The coordinate-\(1\) facets cap both lengths by \(a/|\theta_1|\).  The
other facets give

\[
 \mathbb E\ell_\pm
 \le\min\left({a\over|\theta_1|},
              {2a\over\sum_{j=2}^n|\theta_j|}\right)
 \le {3a\over L}.                                            \tag{5.6}
\]

Therefore

\[
 \mathbb P\{\ell_++\ell_-\le24a/L\}\ge{3\over4}.             \tag{5.7}
\]

The events in (5.5) and (5.7) intersect with probability at least \(1/2\).
On their intersection,

\[
 {\min(\ell_+,\ell_-)\over\ell_++\ell_-}\ge{1\over96}.
\]

Consequently

\[
 {1\over192}\le
 \mathbb E{\min(\ell_+,\ell_-)\over\ell_++\ell_-}
 \le{1\over2}.                                               \tag{5.8}
\]

Equations (5.3) and (5.8) prove

\[
 {\sqrt{12}\over192}J_\theta\le B_\theta\le\sqrt3J_\theta.
                                                                  \tag{5.9}
\]

Averaging and using (5.1) gives the second part of (0.8):

\[
 {c_n\over192}\le\mathbb E_\theta B_\theta\le{c_n\over2}.    \tag{5.10}
\]

The exceptional directions \(\theta_1=0\) have
\(B_\theta=J_\theta=0\), so (5.9) extends to them.

## 6. Three non-extremal stress models

### 6.1 A regular-simplex median cap

Work in

\[
 \Delta_n=\{x_i\ge0:\ \sum_{i=0}^nx_i=1\}
\]

with uniform probability on its affine span.  Let

\[
 E_s=\{x_0\ge s\},\qquad (1-s)^n={1\over2}.                  \tag{6.1}
\]

After the scalar isotropic dilation, its perimeter is

\[
 P(E_s)={n\sqrt n\over2(n+1)\sqrt{n+2}}\,2^{1/n}
       \in\left[{1\over2\sqrt3},{1\over2}\right].            \tag{6.2}
\]

For \(\theta\) in the tangent space \(\{\sum\theta_i=0\}\), a line
through the cap interface has endpoint lengths

\[
 \ell_+=\min_{\theta_i<0}{x_i\over|\theta_i|},\qquad
 \ell_-=\min_{\theta_i>0}{x_i\over|\theta_i|}.               \tag{6.3}
\]

On the interface \(x_0=s\), the remaining coordinates are uniform on a
simplex of mass \(1-s\).  Put

\[
 A_+=\sum_{\substack{i\ge1\\\theta_i>0}}\theta_i,\qquad
 A_-=\sum_{\substack{i\ge1\\\theta_i<0}}|\theta_i|.          \tag{6.4}
\]

If \(\theta_0>0\), then \(A_-=A_++\theta_0\); the other sign is symmetric.
Set \(A=A_-\).  For a subset of random coordinates whose coefficient sum
is \(R\),

\[
 \mathbb P\left\{\min_i{x_i\over|\theta_i|}>t\right\}
   =\left(1-{tR\over1-s}\right)_+^{\,n-1},\qquad
 \mathbb E\min_i{x_i\over|\theta_i|}
   ={1-s\over nR}.                                          \tag{6.5}
\]

The side containing coordinate \(0\) is additionally capped by
\(s/|\theta_0|\).  Since

\[
 \log2\le {ns\over1-s}=n(2^{1/n}-1)\le1,                    \tag{6.6}
\]

the same union-bound argument as for the cube, now at scale
\((1-s)/(nA)\), gives

\[
 \mathbb P\left\{\ell_+,\ell_-\ge{1-s\over8nA}\right\}
 \ge{3\over4},\qquad
 \mathbb E(\ell_++\ell_-)\le{3(1-s)\over nA}.                \tag{6.7}
\]

To verify the second constant, write \(A=a+b\), where
\(b=|\theta_0|\) and \(a\) is the coefficient sum on the same sign side
excluding coordinate \(0\).  The expected length on that side is at most

\[
 {1-s\over n}\min\left({1\over a},
             {{ns/(1-s)}\over b}\right)
 \le {2(1-s)\over n(a+b)},                                  \tag{6.7a}
\]

because either \(a\ge b\) or \(b>a\).  The opposite side has expectation
\((1-s)/(nA)\).  The first constant in (6.7) follows by applying
\(1-(1-z)^{n-1}\le(n-1)z\) separately to the two random coordinate
groups; the deterministic coordinate-\(0\) cap exceeds the displayed
lower scale by (6.6).

Markov's inequality then yields

\[
 \mathbb E{\min(\ell_+,\ell_-)\over\ell_++\ell_-}
 \ge{1\over192}.                                             \tag{6.8}
\]

Because a simplex chord is uniform, (5.3) applies verbatim.  Hence

\[
 {\sqrt{12}\over192}J_\theta\le B_\theta\le\sqrt3J_\theta,
 \qquad
 \mathbb E_\theta B_\theta\asymp n^{-1/2}.                  \tag{6.9}
\]

This cap is a scale test, not an asserted global isoperimetric minimizer.

### 6.2 The isotropic radial-exponential median sphere

Let

\[
 d\mu_n(x)=Z_n^{-1}e^{-\sqrt{n+1}|x|}\,dx.                   \tag{6.10}
\]

Its radial variable has the gamma law
\({\rm Gamma}(n,\sqrt{n+1})\), and the chosen rate makes
\(\operatorname {Cov}(\mu_n)=I\).  Let \(r_n\) be its radial median and
\(E=\{|x|\le r_n\}\).  For \(n\ge2\), the standard gamma-median bounds give

\[
 {n-1/3\over\sqrt{n+1}}\le r_n\le{n\over\sqrt{n+1}}.         \tag{6.11}
\]

By rotation invariance, \(B_\theta\) does not depend on \(\theta\).
Parameterize a boundary point by \(u=N\cdot\theta\).  For \(n\ge32\),
restrict the Crofton integral to

\[
                      {1\over4\sqrt n}\le|u|\le{4\over\sqrt n}.
                                                                  \tag{6.12}
\]
It has a universal positive fraction of the
\(|u|\,d\sigma_{S^{n-1}}\)-mass.  Indeed,
\(c_n\ge(2n)^{-1/2}\), the contribution of
\(|u|<1/(4\sqrt n)\) is at most \(1/(4\sqrt n)\), and
\[
 \mathbb E\{|u|1_{\{|u|>4/\sqrt n\}}\}
 \le {\mathbb Eu^4\over(4/\sqrt n)^3}
 ={3\sqrt n\over64(n+2)}<{1\over20\sqrt n}.                 \tag{6.12a}
\]

On the corresponding lines the
half-length of the spherical chord is

\[
                  r_n|u|\in[c_0,C_0],                       \tag{6.13}
\]

and the conditional line density is proportional to

\[
 t\longmapsto
 \exp\{-\sqrt{n+1}\sqrt{r_n^2(1-u^2)+t^2}\}.                 \tag{6.14}
\]

Equations (6.11)--(6.14), using

\[
 {t^2\over2\sqrt{s^2+t^2}}
 \le\sqrt{s^2+t^2}-s\le{t^2\over2s},                        \tag{6.15}
\]

show with universal constants that its conditional standard deviation is
at most \(C_1\), its normalization is comparable to its density at the
sphere crossings, and both phases have conditional mass at least \(c_2\).
Here are the details needed for uniformity.  On (6.12),
\(r_n|u|\in[c_0,4]\), \(s=r_n\sqrt{1-u^2}\asymp\sqrt n\), and
\(\sqrt{n+1}/s\in[1,3]\).  For \(|t|\le s\), (6.15) sandwiches the
conditional density ratio between Gaussian functions
\(\exp(-3t^2/2)\) and \(\exp(-t^2/(2\sqrt2))\); for \(|t|>s\), its
remaining tail is bounded by a geometric exponential tail.  It follows
that
\[
 c_3q_s(0)\le w_s\le C_3q_s(0),\qquad
 \sigma_s\le C_4.                                           \tag{6.15a}
\]
The interval \([-r_n|u|,r_n|u|]\) and its complement each contain a
fixed subinterval on which \(q_s/q_s(0)\) has a universal positive lower
bound.  Hence their conditional masses are at least \(c_2\).  Also
\(q_s(0)\ge q_s(r_n|u|)\), the density at a sphere crossing.
Consequently

\[
 B_\theta\ge c\,J_\theta,\qquad
 B_\theta\le2\sqrt3\,J_\theta,\qquad
 B_\theta\asymp n^{-1/2}.                                   \tag{6.16}
\]

The radial boundary density is the gamma density at its median.  Stirling's
bounds and (6.11) place it between two positive universal constants, so
\(J_\theta=c_nP(E)\asymp n^{-1/2}\).  The finitely many cases
\(n<32\) obey the same conclusion by the same line estimates with the
fixed range \(|u|\in[1/8,1/2]\) (and \(n=1\) is one-dimensional).

The sphere is a smooth stable high-rank interface, but it is not used as
an exact global minimizer.  Its role is to show that tangential fibers do
not by themselves violate (WFC).

### 6.3 A product-exponential maximum

Let \(Y_1,\ldots,Y_n\) be independent rate-one exponentials; after
translation their product law is isotropic.  Choose \(L\) so that

\[
 E_L=\{\max_iY_i\ge L\},\qquad (1-e^{-L})^n={1\over2}.       \tag{6.17}
\]

For the coordinate direction \(e_j\), a fiber is mixed exactly when all
other coordinates are below \(L\).  Its conditional standard deviation
is one and its minority mass is \(e^{-L}\).  Therefore

\[
 B_{e_j}=J_{e_j}
 =e^{-L}(1-e^{-L})^{n-1}
 ={1\over2}(2^{1/n}-1),                                     \tag{6.18}
\]

and

\[
 P(E_L)=\sum_jJ_{e_j}
 ={n\over2}(2^{1/n}-1)
 \in\left[{\log2\over2},{1\over2}\right].                   \tag{6.19}
\]

This interface has high-rank normal matrix \(Q=I/n\) and perfectly
efficient coordinate fibers.  It is nevertheless not an exact minimizer:
simultaneously beveling the pairwise ridges saves a fixed amount of
perimeter.  Thus it cannot disprove (WFC), but it rules out a proof that
declares high-rank unions of facets extremal solely from their separate
fiber efficiencies.

## 7. The remaining completion-or-saving statement

The fiber argument would close if one proved the following statement
without importing a marginal KLS estimate.

> **Weighted fiber completion-or-saving theorem.**  There are universal
> \(\beta,\eta>0\) with the following property.  Let \(\mu\) be isotropic
> and log-concave and let \(E\) have \(\mu(E)=1/2\).  If
> \[
>             \mathbb E_\theta B_\theta(E)<{\beta\over\sqrt n},
>                                                                  \tag{7.1}
> \]
> then there is a finite-perimeter \(F\) such that
> \[
> \mu(F)=\frac12,\qquad
> P_\mu(F)\le P_\mu(E)-\eta\,\mathbb E_\theta B_\theta(E).
>                                                                  \tag{7.2}
> \]

For an exact minimizer, (7.2) is impossible, so (WFC) follows.  A
quasiminimal version may replace exact mass in (7.2) by

\[
 P_\mu(F)+\psi_\mu\mu(E\mathbin\triangle F)<P_\mu(E).        \tag{7.3}
\]

The proof must resolve all three of the following points simultaneously:

1. choose fiber replacements whose longitudinal saving is quantified by
   \(B_\theta\);
2. bound their internal transverse graph perimeter, not merely the
   cylinder side-wall cost;
3. repair the global volume without assuming an isoperimetric inequality
   for the marginal \(\bar\mu_\theta\).

For compact support there is a fourth point: endpoint-terminated fibers
must either contribute to (7.1) or be assembled into the contact-cell
alternative dictated by (3.16)--(3.17).

The Gaussian and cube calculations show the required scale is sharp.
The cylinder example (3.10) shows that point 2 cannot hold for a single
direction without an averaged rigidity input.

## 8. Circularity audit

1. **KLS-strength of the target.**  Equations (0.4)--(0.6) show explicitly
   that (WFC) proves the dimension-free Cheeger bound.  It may not be
   called a normalization lemma.
2. **Marginal repair.**  The majority base lies in dimension \(n-1\).
   A KLS lower bound for \(\bar\mu_\theta\) has the wrong sign for cheap
   repair, while assuming a balanced base competitor of smaller perimeter
   is the desired saving statement in one lower dimension.  Any induction
   must specify a nondegrading recurrence and an actual construction;
   none is hidden here.
3. **Projected calibration.**  Equation (3.13g) proves that max-flow
   projection yields exactly the quasiminimality budget, with no slack.
   Reusing it as an independent marginal estimate would count the same
   extremality information twice.
4. **Hit-and-run.**  \(I_\theta\) is Boolean conductance for conditional
   resampling.  A dimension-free spectral or conductance estimate for the
   averaged chain would be another KLS formulation and is not assumed.
5. **Directional versus total perimeter.**  The valid inequality is
   \(J_\theta\ge B_\theta/(2\sqrt3)\).  Reversing it requires fiber
   efficiency; integrating (3.6) does not control the transverse term.
6. **Volume repair.**  Invoking a deformation with cost
   \(O(\psi_\mu|\delta|)\) is justified only infinitesimally on a regular
   interface.  A finite deformation with the same bound is precisely the
   quasiminimality comparison (3.2), whose sign does not create a saving.
7. **Test-function completeness.**  No Poincaré inequality, restricted or
   otherwise, is used.
8. **Dimension tracking.**  The only spherical factor is
   \(c_n=\mathbb E|\Theta_1|\le n^{-1/2}\).  All constants
   \(2\sqrt3,192\) are numerical and independent of \(n\).
9. **Degenerate support.**  If \(\mu\) is supported on a \(k\)-dimensional
   affine subspace, all definitions and averages must be taken in that
   subspace with \(k\) replacing \(n\).  Ambient directions orthogonal to
   the support must not be included.  A point mass is excluded.
10. **Support contact.**  Relative-perimeter-free chord endpoints are not
   free in the global replacement.  The exact tensor
   \(\mathcal T_K\) in (3.16) records their aggregate trace.  Treating
   support-killed fibers as a null error would fail already on the cube.

The weighted-fiber route is therefore normalized correctly and survives
the principal model checks, but its completion theorem remains exactly
the new KLS-strength input.
