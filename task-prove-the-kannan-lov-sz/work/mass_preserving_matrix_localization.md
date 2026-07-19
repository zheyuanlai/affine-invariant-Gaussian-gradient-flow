# Mass-preserving matrix localization

## 0. Purpose and status

This note records a controlled stochastic-localization identity which is
different from the usual isotropic-driver bootstrap.  For one fixed Borel set
`S`, the control annihilates the covariance vector of `1_S`.  Consequently
`mu_t(S)` is constant on every path, while all covariance directions except
one are damped.  The exact identities below are proved.  The remaining
terminal long-needle estimate is not proved and is stated explicitly in
Section 4; it is the concurrence/high-filling branch in another language.

Throughout this note the initial measure has a positive smooth log-concave
density with all moments.  The algebraic identities are valid for bounded
predictable controls after the usual stopping.  The feedback rule in Section
2 is discontinuous at a zero of the set-covariance vector and is **not** a
classical SDE there.  A mesh-relaxed weak formulation is required; Section
2.1 records the issue.  Removing the stops and treating nonsmooth or
lower-dimensional measures would still be required in a final proof.

## 1. General controlled localization

Let `p_0` be a probability density on `R^k`.  Let `D_t` be a predictable
positive-semidefinite matrix and let `C_t C_t^T=D_t`.  Define

\[
 p_t(x)={1\over Z_t}
 \exp\{\langle c_t,x\rangle-	frac12\langle B_tx,x\rangle\}p_0(x),
 \qquad dB_t=D_tdt,
\]

where

\[
 dc_t=C_t,dW_t+D_ta_t,dt,
 \qquad a_t=\int x p_t(x)\,dx.
\]

Direct Itô differentiation, including the normalization `Z_t`, gives

\[
                         dp_t(x)
 =p_t(x)\langle x-a_t,C_t,dW_t\rangle.             \tag{1.1}
\]

Put

\[
 A_t=\int(x-a_t)(x-a_t)^T p_t(x)\,dx
\]

and, for a fixed Borel set `S`,

\[
 g_t=\mu_t(S),\qquad
 v_t=\int_S(x-a_t)p_t(x)\,dx.
\]

Integrating (1.1) gives

\[
 dg_t=\langle C_t^Tv_t,dW_t\rangle,
 \qquad d[g]_t=v_t^TD_tv_t\,dt.                     \tag{1.2}
\]

If `T_t(u)=int (x-a_t)(x-a_t)^T <x-a_t,u> p_t(x)dx`,
the covariance equation is

\[
 dA_t=T_t(C_t,dW_t)-A_tD_tA_t\,dt.                 \tag{1.3}
\]

Equations (1.1)--(1.3) use no isoperimetric or spectral-gap input.

There is also a pathwise regression bound.  If `A_t` is nonsingular, then

\[
             v_t^TA_t^{-1}v_t\le g_t(1-g_t).         \tag{1.4}
\]

Indeed, for every vector `q`, Cauchy--Schwarz gives
`<q,v_t>^2 <= g_t(1-g_t) q^TA_tq`; optimize in `q`.

## 2. The mass-preserving control

Assume `A_t` is positive definite.  If `v_t!=0`, set

\[
 w_t=A_t^{-1/2}v_t,
 \quad P_t=I-{w_tw_t^T\over |w_t|^2},
 \quad D_t=A_t^{-1/2}P_tA_t^{-1/2}.                 \tag{2.1}
\]

At `v_t=0`, choose provisionally any rank-`k-1` orthogonal projection
`P_t` (and take `P_t=0` when `k=1`).  Then

\[
                         D_tv_t=0.                   \tag{2.2}
\]

Thus (1.2) proves the exact pathwise conservation law

\[
                         g_t=g_0.                    \tag{2.3}
\]

For `v_t!=0`, the covariance drift in (1.3) is

\[
 -A_tD_tA_t
 =-A_t+{v_tv_t^T\over v_t^TA_t^{-1}v_t}.            \tag{2.4}
\]

It has rank `k-1` in covariance-normalized coordinates.  In particular,

\[
 \operatorname {tr}(A_tD_t)=k-1,
 \qquad
 {d\over dt}\big|_{\rm drift}\log\det A_t=-(k-1). \tag{2.5}
\]

The second identity follows from
`tr(A_t^{-1}A_tD_tA_t)=tr(D_tA_t)=k-1`; the Itô correction
to `log det A_t` is nonpositive.  Hence this control contracts covariance
volume while spending no quadratic variation of the chosen set mass.

For comparison, the full whitening control `D_t=A_t^{-1}` has covariance
drift `-A_t` but spends at most `g_t(1-g_t)dt` of set-mass quadratic
variation by (1.4).  Formula (2.1) deletes precisely the one whitened mode
which causes that expenditure.

### 2.1 Zero-set well-posedness and the relaxed control

The pointwise convention above does not by itself define a strong or weak
feedback SDE.  The tempting convention `P=I` at `v=0` is actually
inconsistent.  In dimension one, take `mu_0=N(0,1)` and a symmetric interval
`S=[-r,r]`.  For a Gaussian posterior `N(a,A)`,

\[
 v=A\,\partial_a\mu_{a,A}(S),
\]

so `v=0` exactly when `a=0`.  The rule `D=A^{-1}` at `v=0` and `D=0`
elsewhere would give a continuous natural parameter `c` with

\[
 d[c]_t=A_t^{-1}1_{\{c_t=0\}}dt.
\]

The occupation-density formula gives
`int 1_{c_t=0}d[c]_t=0`, hence `[c]=0` and `c=0`; the displayed SDE then
forces strictly positive quadratic variation.  Thus no continuous weak
solution exists for that convention.

A legitimate formulation follows the relaxed construction for
mass-preserving rank-one localization.  On a time mesh, freeze during each
step a rank-`k-1` whitened projection annihilating the preceding value of
`w`; at `w=0` choose any such projection.  After covariance and parameter
stopping the laws are tight.  Subsequential limits give a predictable
positive contraction `P_t` satisfying

\[
 0\preceq P_t\preceq I,
 \qquad \operatorname {tr}P_t=k-1,
 \qquad P_tA_t^{-1/2}v_t=0                         \tag{2.6}
\]

almost everywhere, and
`D_t=A_t^{-1/2}P_tA_t^{-1/2}`.  Equations (2.2), (2.3), and the trace
identity in (2.5) remain valid in integrated form; (2.4) applies whenever
`v_t!=0`.  This requires a full tightness proof on compactly
truncated measures and a later truncation limit.  More importantly, a
terminal estimate must hold uniformly over every such relaxed limit; no
canonical pathwise control has yet been proved to exist.

## 3. Perimeter remains a nonnegative local martingale

Suppose temporarily that `S` has a compact `C^1` boundary.  Its weighted
perimeter under `p_t` is

\[
 P_t(S)=\int_{\partial S}p_t\,d\mathcal H^{k-1}.
\]

Integrating (1.1) over the fixed boundary yields

\[
 dP_t(S)=\left\langle C_t^T
 \int_{\partial S}(x-a_t)p_t(x)\,d\mathcal H^{k-1}(x),dW_t\right\rangle.
                                                               \tag{3.1}
\]

Thus stopped perimeter is a martingale and unrestricted perimeter is a
nonnegative local martingale, so

\[
                         E P_t(S)\le P_0(S).          \tag{3.2}
\]

Approximation by smooth sets gives the same one-sided statement for finite
perimeter sets at every bounded stopping time for which the surface
integrals are uniformly integrable.

If a relaxed process converges to a one-dimensional log-concave needle
`nu_omega` on a line of direction `u_omega`, then (2.3) says

\[
                       \nu_\omega(S)=g_0             \tag{3.3}
\]

for almost every terminal path.  The collapse of transverse directions
does not make ambient perimeter smaller than the one-dimensional Minkowski
boundary of the trace.  For a transverse smooth crossing this is the
identity

\[
 \lim P_t(S)=\sum_{x\in\partial S\cap L_\omega}
 {\rho_\omega(x)\over |n_S(x)\cdot u_\omega|}
 \ge \nu_\omega^+(S\cap L_\omega),                 \tag{3.4}
\]

and tangencies are obtained by lower semicontinuity.  Therefore the
one-dimensional log-concave Cheeger inequality would give

\[
 P_0(S)\ge c\,g_0(1-g_0)
 E\,{1\over \sigma_\omega},                         \tag{3.5}
\]

where `sigma_omega^2` is the variance of the terminal needle.  This last
passage requires a full convergence and uniform-integrability proof; (3.5)
is recorded only as the precise intended terminal reduction.

## 4. The exact remaining obstruction

Disintegration of the initial covariance gives

\[
 \operatorname {Cov}(\mu_0)
 =\operatorname {Cov}(b_\omega)
   +E\,[\sigma_\omega^2u_\omega u_\omega^T].         \tag{4.1}
\]

For isotropic `mu_0`, (4.1) implies the operator inequality

\[
              E[\sigma_\omega^2u_\omega u_\omega^T]\preceq I. \tag{4.2}
\]

It does **not** algebraically imply `E sigma_omega^2<=C`: taking
`sigma_omega^2=k` and an isotropic uniform direction saturates (4.2).
Consequently neither (4.2) nor Jensen closes (3.5).

The missing statement can be isolated as follows.

> **Terminal-needle lemma (unproved).**  For every mesh-relaxed limit of
> (2.1), started from an isotropic log-concave probability and any Borel
> `S` with `0<mu(S)<1`, the terminal needles satisfy
> `E sigma_omega^2<=C`, or at least
> `E(1/sigma_omega)>=c`, with universal constants.

This lemma is not a consequence of covariance disintegration alone.  A
putative bad family must have long needles with directions of large
effective rank and with almost coincident barycenters, because (4.1) leaves
little barycenter covariance when the second term nearly saturates `I`.
That is exactly the concurrent/orthogonal-radial branch exposed by the
large-interior-filling analysis.  For a Euclidean ball, central chords have
variance of order `k`, but the Gaussian-conditioning disintegration chooses
typical chords near the boundary and their conditional variance is of order
one; the rare central chord cannot be used as a pointwise counterexample.

The following tempting shortcuts are invalid.

1. Taking the trace of (4.2) gives only `E sigma^2<=k`.
2. Selecting a terminal direction after seeing the localization path is
   not bounded by the initial covariance operator; it is an adaptive
   selection problem.
3. Pointwise conditional variances of a log-concave measure can exceed the
   corresponding marginal variance, so no elementary conditioning
   monotonicity proves the terminal-needle lemma.
4. Claiming that `B_t` becomes uniformly positive contradicts (2.3): if all
   directions collapsed to a point, the nontrivial fixed mass `g_0` could
   not remain pathwise constant.  A surviving needle is necessary.

Thus the control removes the standard stochastic-localization covariance-
maximum obstruction but replaces it with one explicit geometric inverse.
Proving the terminal-needle lemma, or proving a thin-shell estimate for its
high-rank concurrent disintegration, would close this route; assuming it is
conjecture-strength.
