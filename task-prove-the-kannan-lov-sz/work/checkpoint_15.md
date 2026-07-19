# Checkpoint 15: convolution amplification and the weighted-slice gate

## 1. Candidate-proof status

There is still no complete proof of KLS, and no KLS conclusion is claimed
at this checkpoint.  The new work proves a dimension-free forward
self-convolution theorem, closes the weighted-Fisher estimate for every
pure translating slice family, and proves a full weighted-Fisher estimate
for smooth full-support Gaussian conditional fibers after unit transverse
Gaussian noise.  It identifies two exact global gates that survive the
present audits and supplies explicit counterexamples to three tempting
reverse or slice-local arguments.

Throughout this checkpoint, every covariance normalization is Euclidean.
Smooth calculations are asserted for the regularity class stated in the
underlying notes; hard-support passages are claimed only where a direct
finite-difference or closed-form argument is supplied.

## 2. Strict spectral amplification under normalized self-convolution

Let \(X,Y\) be independent with the same centered isotropic law \(\mu\),
put

\[
 S=\frac{X+Y}{\sqrt2},\qquad \nu=\mathcal L(S),
\]

and write

\[
 a=\lambda _1(\mu)=C_P(\mu)^{-1},\qquad
 b=\lambda _1(\nu)=C_P(\nu)^{-1}.
\]

The proof in `self_convolution_gradient_rigidity.md` does not assume that
the bottom of either spectrum is attained.  It gives

\[
 \boxed{\quad b(1-b)\le4(b-a),\qquad
 b\ge\frac{\sqrt{9+16a}-3}{2}.\quad}                 \tag{2.1}
\]

The algebra behind (2.1) is useful.  For a centered unit test function
\(f(S)\), let

\[
 h(X)=E[f(S)\mid X],\qquad
 R=f(S)-h(X)-h(Y).
\]

If \(q=E|\nabla f|^2\), then the product Poincare inequality on the
two-coordinate-centered residual gives

\[
 E_R:=E|\nabla_{X,Y}R|^2\le2(q-a).                    \tag{2.2}
\]

For \(Z=\nabla f(S)\), \(U=E[Z\mid X]\), \(V=E[Z\mid Y]\), and
\(W=Z-U-V+EZ\), the vector Hoeffding decomposition gives

\[
 \operatorname {Var}(Z)=2\operatorname {Var}(U)+E|W|^2,
 \qquad E_R=\operatorname {Var}(U)+E|W|^2,            \tag{2.3}
\]

and hence

\[
 \operatorname {Var}(\nabla f)\le2E_R\le4(q-a).      \tag{2.4}
\]

Choosing \(f\) in a spectral window \([b,b+\varepsilon]\), testing the
weak generator equation against coordinates, and sending
\(\varepsilon\downarrow0\) gives
\(\operatorname {Var}(\nabla f)\ge b(1-b)-o(1)\), which proves (2.1)
even at a continuous spectral edge.

If \(a_k\) denotes the gap after \(k\) normalized dyadic convolutions,
then \(a_k\uparrow1\).  When \(a_{k+1}\le1/2\),

\[
 a_{k+1}\ge\frac87a_k.                               \tag{2.5}
\]

This is only a forward amplification theorem.  Reaching a fixed gap can
take \(O(\log(1/a_0))\) steps, so (2.1) supplies no lower bound on the
initial gap.

## 3. Reverse mechanisms that fail, and the saturation structure

### 3.1 Conditional projection is not an energy contraction

Let \(f\) be an exact normalized bottom mode of \(\mu\), with eigenvalue
\(a\), and set

\[
 H=\frac{f(X)+f(Y)}{\sqrt2},\qquad
 g=E[H\mid S],\qquad r=H-g.
\]

Then

\[
 \boxed{\quad
 \mathcal E_\nu(g)=a(1-2\|r\|_2^2)+\mathcal E_{\mu^2}(r).
 \quad}                                               \tag{3.1}
\]

There is no general bound \(\mathcal E_\nu(g)\le\mathcal E_{\mu^2}(H)\).
For \(\mu=\operatorname {Unif}[-1,1]\), take

\[
 f(x)=\sqrt2\sin(\pi x/2),\qquad a=\pi^2/4.
\]

Writing \(t=\sqrt2s\), direct conditional averaging gives

\[
 g(s)=\frac4\pi\frac{\sin(\pi t/2)}{2-|t|},\qquad |t|<2.
\]

With \(J=\int_0^\pi\sin^2(v)\,dv/v\),

\[
 \operatorname {Var}_\nu(g)=\frac{8J}{\pi^2},\qquad
 \mathcal E_\nu(g)=4(J-1/2)>\frac{\pi^2}{4}.          \tag{3.2}
\]

The strict inequality follows from the elementary lower bound

\[
 J\ge\frac{\pi^2}{8}-\frac{\pi^4}{192}
 +\frac{\pi^6}{13824}+\frac5{16}+\frac1{4\pi^2}
 >\frac12+\frac{\pi^2}{16}.                           \tag{3.3}
\]

Thus discarding the sum-fiber residual cannot reverse (2.1).

### 3.2 Scalar convolution Lyapunovs do not supply the reverse budget

For the standardized Gamma law
\(X_k=(\Gamma(k,1)-k)/\sqrt{k}\), normalized self-convolution sends
\(X_k\) to \(X_{2k}\), while its exact Euclidean gap is

\[
 g_k=\frac{k^2}{(k+1)^2}.                              \tag{3.4}
\]

This calibration exposes the following obstructions.

* Raw Fisher information equals \(J(X_k)=k/(k-2)\) for \(k>2\), but is
  infinite at the mandatory hard-support limits \(k\le2\).
* Gaussian regularizations of the one-sided exponential retain the gap
  change \(1/4\to4/9\), while every bounded inverse-Fisher decrement
  tends to zero.
* For each fixed even \(m\ge4\), the strongly log-concave Hermite
  perturbation
  \(p_{\varepsilon,m}\propto
  \exp[-x^2/2-\varepsilon H_m(x)/\sqrt{m!}]\), after isotropization,
  satisfies
  \[
  a=1-\frac{m}{m-2}\varepsilon^2+O_m(\varepsilon^3),\qquad
  b=1-\frac{m}{m-2}2^{2-m}\varepsilon^2
       +O_m(\varepsilon^3).
  \]
  For the unit-noise Fisher/MMSE deficit \(\mathcal F\),
  \[
  \lim_{\varepsilon\downarrow0}
  \frac{\log(b/a)}{\mathcal F(\mu)-\mathcal F(\mathsf T\mu)}
  =\frac{2^m}{m-2}.                                    \tag{3.5}
  \]
  Hence no fixed-noise scalar Fisher charge is universal.
* Entropy deficit and Hadamard sum/difference mutual information survive
  these one-dimensional tests, but they are additive under products.
  The raw quantities fail a universal budget on tensor powers, division
  by dimension fails under Gaussian padding, and bounded transforms
  saturate.

These examples rule out the tested scalar budgets; they do not rule out a
directional, operator-valued, or function-specific Lyapunov.

### 3.3 A coherent mean Hessian is forced at forward saturation

In the smooth convex setting, suppose an attained exact bottom mode of
\(\nu\) has eigenvalue \(b\).  With the notation in (2.3), put

\[
 v=\operatorname {Var}(U),\quad w=E|W|^2,\quad
 B=E\|D^2f\|_{\mathrm {HS}}^2,\quad C=E D^2f.
\]

Then

\[
 \operatorname {Var}(\nabla f)=2v+w,\quad E_R=v+w,
 \quad B\le b^2,
\]

\[
 4av\le B+\|C\|_{\mathrm {HS}}^2,qquad B\ge2aE_R.   \tag{3.6}
\]

Consequently, along any such sequence with \(a\to0\) and
\(b/a\to4/3\),

\[
 \liminf
 \frac{\|E D^2f\|_{\mathrm {HS}}^2}
      {E\|D^2f\|_{\mathrm {HS}}^2}\ge\frac12.        \tag{3.7}
\]

The spectral-window and approximation versions of (3.7) have not been
proved.  In its present scope, (3.7) identifies a coherent mean-Hessian,
equivalently directional third-moment, block that any reverse rigidity
argument must exclude.

## 4. The Gaussian-channel affine gate

Let \(X\sim\mu\) be centered and isotropic, let \(Y=X+G\) with
\(G\sim N(0,I)\), and let

\[
 h(y)=E[g(X)\mid Y=y],\qquad g\perp\{1,X_1,\ldots,X_n\}.
\]

Affine orthogonality and the output score give the unconditional bound

\[
 |E\nabla h|\le2^{-1/2}\|g\|_2.                       \tag{4.1}
\]

For the posterior law \(\mu_y\), define

\[
 A_y=\operatorname {Cov}_y(X),\qquad
 a_y=E_y\nabla g,
\]

and subtract from \(g\) its posterior mean and affine predictor to obtain
\(r_y\).  Posterior strong convexity gives the exact decomposition

\[
 E_\mu\nabla g-E_q\nabla h
 =E[(I-A_Y)a_Y]-E\operatorname {Cov}_Y(r_Y,X),          \tag{4.2}
\]

with

\[
 \left|E\operatorname {Cov}_Y(r_Y,X)\right|
 \le\|D^2g\|_{L^2(\mu)}.                              \tag{4.3}
\]

All non-affine posterior modes are therefore controlled.  The exact
remaining term is

\[
 \boxed{\quad
 \left|E[(I-A_Y)a_Y]\right|
 \le C\|D^2g\|_{L^2(\mu)}.
 \quad}                                               \tag{AFF}
\]

A pointwise posterior version is false because a posterior-affine
function has zero Hessian.  The global estimate (AFF) remains unproved and
is KLS-strength; invoking a dimension-free Poincare inequality for the
output law would be circular.

## 5. Exact weighted-Fisher and reinforced-Prekopa ledger

Let

\[
 r(s,z)=e^{-U(s,z)},\quad \rho(s)=\int r(s,z)\,dz=e^{-\Phi(s)},
 \quad q_s=r/\rho,
\]

and, in the regular setting, solve

\[
 L_sg_s=\ell_s:=\partial_s\log q_s,qquad
 F_s=\nabla_zg_s,qquad H_s=D^2_{zz}U.
\]

Put

\[
 C_s=E_s\|D_zF_s\|_{\mathrm {HS}}^2,qquad
 B_s=E_s\langle H_sF_s,F_s\rangle.
\]

Conditional Bochner gives

\[
 I_\perp(s):=E_s\ell_s^2=C_s+B_s.                    \tag{5.1}
\]

With

\[
 \alpha_s=\partial_sg_s-\tfrac12|F_s|^2,qquad
 \mathcal Q_s=D^2U[(1,-F_s),(1,-F_s)],
\]

direct differentiation gives

\[
 L_s\alpha_s=\Phi''-\mathcal Q_s-
                  \|D_zF_s\|_{\mathrm {HS}}^2,
\qquad
 \Phi''=E_s\mathcal Q_s+C_s.                         \tag{5.2}
\]

If \(\tau\) is the canonical Stein kernel of the centered variance-one
marginal \(\rho\), then

\[
 \int\rho\tau^2(C_s+E_s\mathcal Q_s)\,ds
 =\int\rho\tau^2\Phi''\,ds
 =1-E(\tau')^2\le1.                                  \tag{5.3}
\]

Thus \(C_s\) and the material curvature have a universal weighted budget.
The unresolved term is

\[
 \mathcal B=\int\rho(s)\tau(s)^2B_s\,ds.             \tag{5.4}
\]

The full-space Hodge identity reproduces exactly (5.3): the apparent
\(B_s\) contribution in the Hessian quadratic form cancels against the
horizontal and mixed terms.  Reinforced Prekopa alone therefore does not
bound (5.4).

The material identities also give, for the conditional centroid
\(m(s)=E_sZ\),

\[
 m'=-E_sF_s,qquad
 m''=-E_s[(Z-m)(\mathcal Q_s+\|D_zF_s\|_{\mathrm {HS}}^2)],
\qquad E[\tau(S)m'(S)]=0.                             \tag{5.5}
\]

These are the global-in-\(s\) data that a valid proof of (5.4) must use.

## 6. Complete positive class: every pure translating slice family

Suppose the unsmoothed conditional laws have the form

\[
 p_s(x)=Z_\varphi^{-1}e^{-\varphi(x-m(s))}
\]

and the joint law is isotropic and log-concave.  If \(Y\) has density
proportional to \(e^{-\varphi}\), then
\(\operatorname {Cov}(Y)\preceq I\).  For

\[
 K=\overline{\operatorname {conv}}\{\nabla\varphi(y)\},
\]

with all subgradients and normal cones included in the nonsmooth case,
one has the dimension-free inradius bound

\[
 K\supset\frac3{16}B_2^d.                             \tag{6.1}
\]

Indeed, every one-dimensional marginal has variance at most one, hence
density maximum at least \(3/16\); its extreme marginal score is at least
that maximum.  Testing joint convexity on \((1,m'(s))\) now gives, in
distributions,

\[
 |Dm'|\le\frac{16}{3}D\Phi'.                          \tag{6.2}
\]

At a mode, write \(b=m'\).  Cross-isotropy, the one-dimensional Stein
identities, \(E\tau^2\le400\), and
\(E[\tau^2(\Phi')^2]\le4\) yield

\[
 E[\tau^2|m'|^2]
 \le3208\left(\frac{16}{3}\right)^2.                 \tag{6.3}
\]

After unit transverse Gaussian smoothing, the continuity velocity is
\(m'(s)\), and conditional Jensen gives

\[
 I_\perp(s)\le|m'(s)|^2.                              \tag{6.4}
\]

Equations (6.3)--(6.4) prove the desired weighted-Fisher inequality for
every pure translating family, uniformly in the transverse dimension.
The direct finite-difference proof includes asymmetric exponential fibers
and hard endpoints.  This theorem depends on a score-range **inradius**;
a score-radius upper bound would have the wrong direction.

## 7. Pointwise post-Gaussian curvature--Korn is false

For \(q_a=\operatorname {Unif}[-a,a]*\gamma\) and \(F_a(x)=x/a\),

\[
 E_{q_a}|F_a'|^2=a^{-2},\qquad
 E_{q_a}[U_a''F_a^2]\ge\frac1{64a}.                  \tag{7.1}
\]

More strongly, the same boundary-layer obstruction occurs for the exact
centered Poisson field of one fixed isotropic log-concave law.  Start from

\[
 p(r,x)=\tfrac12e^{-r}\mathbf1_{\{r>0,\ |x|<r\}},
\]

put \(S=(R-2)/\sqrt2\), \(Y=X/\sqrt2\), and
\(a=S+\sqrt2\).  Then \(Y\mid S=s\) is uniform on \([-a,a]\).  After unit
Gaussian smoothing, its exact Poisson field satisfies

\[
 E_sF_s=0,qquad C_s\le a^{-2},\qquad B_s\ge(256a)^{-1}.
                                                               \tag{7.2}
\]

Thus no universal slice-by-slice estimate
\(B_s\lesssim C_s+|E_sF_s|^2\) exists.  The integrated target survives:
the scalar Stein kernel is \(\tau(S)=a/\sqrt2\), \(B_s\le1\), and

\[
 \int\rho\tau^2B_s\,ds\le\frac32.                   \tag{7.3}
\]

This example forces any proof to retain the joint scalar geometry and
Stein weighting.

## 8. Shape-changing centroid and Gaussian-fiber WFI results

Two shape-changing calculations are verified.  The first now gives a full
weighted-Fisher theorem for a smooth post-noise subclass; the second remains
a geometric hard-support result.

1. For a smooth full-support conditional Gaussian family
   \[
   U(s,z)=W(s)+\tfrac12(z-m(s))^TQ(s)(z-m(s)),\qquad R=Q^{-1},
   \]
   joint convexity implies \(R''\preceq0\) and the Schur condition
   \[
   W''-m''\cdot x+\tfrac12x^T(-R'')x\ge0.
   \]
   If \(S\) is centered and variance one,
   \(E[S\,m(S)]=0\) as an absolutely defined vector moment, and
   \(ER\preceq\Lambda I\), then
   \[
   E[\tau(S)^2|m'(S)|^2]\le C_{\rm cent}\Lambda.    \tag{8.1}
   \]
   The proof uses the corrected one-dimensional hazard estimates
   \[
   \tau(s)\le400(1+|s|),\qquad
   \tau(s)|\Phi'(s)|\le400(1+|s|).
   \]
   The earlier constant-one pointwise score claim is not verified and is
   not used.

   If in addition \(R(s)\succeq I\), as supplied by unit transverse
   Gaussian noise, the exact conditional Poisson field is
   \[
   F_s=-m'(s)+A_s(z-m(s)),\qquad A_sR+RA_s=-R'.
   \]
   With
   \[
   C_s=\|A_s\|_{\rm HS}^2,\qquad
   B_s=m'^TR^{-1}m'+B_s^0,\qquad
   B_s^0=\operatorname {tr}(R^{-1}A_sRA_s),
   \]
   diagonalizing \(R(s)\) at a fixed \(s\) and pairing the two ordered
   off-diagonal terms gives the exact identity
   \[
   \boxed{\quad C_s+B_s^0=J_s,
   \qquad
   J_s:=\tfrac12\operatorname {tr}
   (R^{-1}R'R^{-1}R')\le\Phi''(s).\quad}              \tag{8.2}
   \]
   Therefore the scalar Stein-curvature identity yields
   \[
   \int\rho\tau^2(C_s+B_s^0)\,ds\le1.
   \]
   Since \(R\succeq I\) gives \(R^{-1}\preceq I\), (8.1) controls the
   centroid term and proves the dimension-free smooth Gaussian-fiber bound
   \[
   \boxed{\quad
   \int\rho(s)\tau(s)^2(C_s+B_s)\,ds
   \le1+C_{\rm cent}\Lambda.\quad}                   \tag{8.3}
   \]
   This calculation allows noncommuting, rotating covariance matrices.
   It is a \(C^2\), full-support, unit-noise subclass theorem only; no
   nonsmooth, distributional, or zero-noise passage is claimed.  The lower
   bound \(R\succeq I\) is essential at the last comparison; more generally
   \(R\succeq\kappa I\) costs \(C_{\rm cent}\Lambda/\kappa\).
2. For planar uniform interval fibers
   \(l(s)\le z\le u(s)\), convexity gives a dimension-free weighted
   bound on the centroid velocity, including wedge curvature atoms.  This
   is a hard-support geometric statement; its identification with the
   post-noise Poisson field has not been proved.

For a general regular slice, set \(b=m'(s)\) and
\(F_s^0=F_s+b\), so \(E_sF_s^0=0\).  A fixed affine shear and conditional
Bochner identity give the exact formula

\[
 \boxed{\quad
 E_sD^2U[(1,m'(s)),(1,m'(s))]
 =\Phi''(s)+C_s+B_s^0,
 \quad}                                               \tag{8.4}
\]

where
\(B_s^0=E_s\langle H_sF_s^0,F_s^0\rangle\).  Thus the affine-shear
curvature does not remove the nonlinear residual; it reproduces it
exactly.

## 9. Exact surviving gates

The weighted route would close if one proved, with a universal constant,

\[
 \int\rho\tau^2B_s
 \le C\int\rho\tau^2(C_s+E_s\mathcal Q_s).            \tag{WK}
\]

The cone example shows that (WK) must be integrated in \(s\); the
translation theorem shows that all pure centroid motion is already under
control.  The unresolved sector is genuinely nonlinear shape change.

Independently, the near-linear Gaussian-channel route requires (AFF).
Even a proof of (WK), WFI, or the corresponding directional-MMSE estimate
would presently close only the near-linear spectral branch.  The
tensor-extremizer program still lacks a proved dimension-free dichotomy
forcing every hypothetical small-gap mode into that branch.  Accordingly,
the two load-bearing tasks are:

1. prove the integrated residual estimate (WK), without a conditional
   Poincare insertion;
2. prove the nonlinear/near-linear dichotomy, or find a different argument
   that applies the weighted estimate to every low mode.

Neither statement is assumed in this checkpoint.
