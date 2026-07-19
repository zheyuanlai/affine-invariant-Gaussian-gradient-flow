# Mass-preserving localization and the rank-one curvature reduction

This note records a dimension-free part of a set-specific localization route.
It deliberately stops at the one adaptive posterior variance which is not yet
controlled.  In particular, no operator-norm estimate for a localization
covariance is assumed.

## 1. An all-but-one-direction curvature lemma

Let `mu` be a log-concave probability on `R^n` with density `exp(-V)`, and
let `u` be a unit vector.  Assume, in the distributional sense, that

\[
 V(x)-{\alpha\over2}|P_{u^\perp}x|^2
\quad\hbox{is convex}                                      \tag{1}
\]

for some `alpha>0`.  Put

\[
 s^2=\operatorname {Var}_\mu\langle X,u\rangle.
\]

Then there is a universal `c>0` such that

\[
 \boxed{\quad \psi_\mu\ge
 {c\over\sqrt{\alpha^{-1}+s^2}}.\quad}                    \tag{2}
\]

Here is a needle proof, conditional only on the standard balanced Euclidean
needle-disintegration theorem.  Fix a finite-perimeter set `A`, write
`p=mu(A)`, and disintegrate `mu` into one-dimensional log-concave needles
`mu_omega`, supported on affine lines with unit directions `theta_omega`, so
that

\[
 \mu_\omega(A)=p                                             \tag{3}
\]

for almost every `omega`.  The disintegration preserves integrals, and
Minkowski enlargement along each line plus Fatou gives

\[
 \mu^+(A)\ge\int\mu_\omega^+(A)d\pi(\omega).               \tag{4}
\]

The affine factor introduced by Euclidean localization is log-affine (or has
additional log-concavity), so (1) implies that the one-dimensional potential
on the needle has second derivative at least

\[
 \beta_\omega=\alpha(1-a_\omega^2),\qquad
 a_\omega=|\langle\theta_\omega,u\rangle|.                 \tag{5}
\]

Let `tau_omega^2` be the variance of the affine line coordinate and put

\[
 r_\omega^2=\operatorname {Var}_{\mu_\omega}
               \langle X,u\rangle
             =a_\omega^2\tau_\omega^2.                    \tag{6}
\]

The one-dimensional log-concave Cheeger estimate and the one-dimensional
strong-convexity estimate give, with universal constants,

\[
 \psi_{\mu_\omega}\ge c_0
 \max\left({1\over\tau_\omega},
            \sqrt\alpha\sqrt{1-a_\omega^2}\right).        \tag{7}
\]

If `a_omega^2>=1/2`, the first term in (7) is
`a_omega/r_omega>=1/(sqrt(2)r_omega)`.  If
`a_omega^2<1/2`, the second is at least `sqrt(alpha/2)`.
In either case,

\[
 \psi_{\mu_\omega}\ge
 {c_0/\sqrt2\over\sqrt{\alpha^{-1}+r_\omega^2}}.          \tag{8}
\]

Equations (3), (4), and (8) imply

\[
 {\mu^+(A)\over\min(p,1-p)}
 \ge c_1\int {d\pi(\omega)\over
                  \sqrt{\alpha^{-1}+r_\omega^2}}.         \tag{9}
\]

The function `x -> (alpha^{-1}+x)^{-1/2}` is convex.  Moreover, the
conditional-variance identity gives

\[
 \int r_\omega^2d\pi(\omega)
 \le\operatorname {Var}_\mu\langle X,u\rangle=s^2.        \tag{10}
\]

Jensen's inequality in (9) proves (2).  The same argument works on the
minimal affine support.  Approximation of the density and of finite-perimeter
sets is still to be written if (2) is used in a final proof; no constant in
the displayed argument depends on the dimension.

## 2. A set-mass-preserving localization

Let `S` be a fixed Borel set of mass `p` and let `p_t` be the stochastic
localization posterior

\[
 p_t(x)={1\over Z_t}\exp\left(c_t\cdot x-
                    {1\over2}x^TQ_tx\right)p_0(x).         \tag{11}
\]

With `m_t=E_{p_t}X`, put

\[
 v_t=\operatorname {Cov}_{p_t}(1_S,X),\qquad
 C_t=P_{v_t^\perp},                                       \tag{12}
\]

where `C_t=I` when `v_t=0`.  In the standard martingale
parametrization,

\[
 dp_t(x)=p_t(x)(x-m_t)^TC_t\,dW_t,
 \qquad dQ_t=C_t^2dt.                                     \tag{13}
\]

Consequently

\[
 dp_t(S)=v_t^TC_t\,dW_t=0,                                \tag{14}
\]

so `p_t(S)=p` pathwise.  If `S` has a sufficiently regular fixed boundary,
pointwise martingality and Tonelli give

\[
 \mathbb E\,p_T^+(S)=p_0^+(S).                            \tag{15}
\]

For a general finite-perimeter set, (15) requires the usual weighted-BV
approximation and uniform-integrability argument; equality can be weakened
to the direction needed below.

Writing `e_t=v_t/|v_t|` whenever `v_t` is nonzero, one has

\[
 Q_T=TI-\int_0^T e_te_t^Tdt                               \tag{16}
\]

apart from intervals on which `v_t=0`, which only add full curvature.

## 3. The deterministic rank-one spectral fact

Set

\[
 R_T=\int_0^T e_te_t^Tdt.
\]

Then `R_T` is positive semidefinite and `Tr R_T<=T`.  If its eigenvalues are
`rho_1>=rho_2>=...`, then

\[
 \rho_2\le \sum_{j\ge2}\rho_j
          \le T-\rho_1.                                   \tag{17}
\]

Since the eigenvalues of `Q_T=TI-R_T` in the same eigenbasis are
`q_j=T-rho_j`, (17) shows

\[
 \boxed{\quad q_2(Q_T)\ge T/2.\quad}                      \tag{18}
\]

Indeed, if `q_1=T-rho_1<=T/2`, then
`q_2>=T-q_1>=T/2`; if `q_1>T/2`, all eigenvalues exceed
`T/2`.  Thus, for every sample path, there is a unit eigenvector `u_T` such
that

\[
 Q_T\succeq {T\over2}P_{u_T^\perp}.                       \tag{19}
\]

This conclusion is unaffected by arbitrary rotation of the protected
directions `e_t`.

## 4. The exact remaining random quantity

The endpoint potential is convex plus `x^TQ_Tx/2`; hence (2), (14), and
(19) give

\[
 p_T^+(S)\ge c\min(p,1-p)
 \left({2\over T}+
 \operatorname {Var}_{p_T}\langle X,u_T\rangle\right)^{-1/2}.
                                                                    \tag{20}
\]

Averaging and using (15) yields the rigorous reduction

\[
 \boxed{\quad p_0^+(S)\ge c\min(p,1-p)\,
 \mathbb E\left[
 {1\over\sqrt{2/T+\operatorname {Var}_{p_T}
                         \langle X,u_T\rangle}}
 \right].\quad}                                           \tag{21}
\]

Therefore this route would prove KLS if, for one universal `T>0`, one could
show

\[
 \mathbb P\left\{
 \operatorname {Var}_{p_T}\langle X,u_T\rangle\le C
 \right\}\ge c                                             \tag{22}
\]

for every isotropic log-concave initial law and every set `S`.

The eigenvector `u_T` in (22) is path-dependent.  The ordinary law of total
covariance controls `E u^TA_Tu` for each deterministic `u`, but does not
control `E u_T^TA_Tu_T`; replacing the latter by `E||A_T||_op` is exactly the
known covariance-process obstruction.  Thus (22) is not presently proved by
the reduction above.  Any continuation must exploit that `u_T` arises from
the protected binary correlation `v_t`, rather than treating it as an
arbitrary adaptively selected direction.

