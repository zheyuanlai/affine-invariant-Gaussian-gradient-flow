# Clean-room audit of `fixed_scale_wedge_extremality.md`

## Verdict

The numerical chain in Sections 1--4 is correct except for a scale error in
(3.16) and two zero-set/endpoint conventions.  In particular, the conversion
from the endpoint estimate for `A_*` to `A_*/J`, the central-good Fisher trace,
the effective-rank seed, the derivative identity (4.4), and the powers in the
wedge comparison all check out.  The Laplace example in Section 5.1 is valid.

The assertion in Section 6 that the implication is already completely
rigorous needs three small logical repairs.  More substantively, the proposed
phase-charge lemma is stated uniformly in a trace threshold `r` which never
appears in its conclusion.  That formulation is substantially stronger than
the version actually used and is vulnerable to an arbitrarily small,
high-rank phase packet added to an affine configuration.

## 1. Verified constants in Sections 1--2

1. In (1.8),
   `c_G/8=I_0/4` and `I_0/4<1/(16I_0)` because
   `4I_0^2=2/pi<1`.  Hence the displayed value `c_0=I_0/4` is correct.

2. From `K<=4/psi^2` and `p(1-beta)<=psi/2`, (2.4) is correct.  Substitution
   of `s=alpha K` into the resampling lower bound and the profile upper bound
   gives every constant in (2.5).

3. Monotonicity of `sqrt(r)A_*(r)` gives

   \[
   H(s)-H(2s)\ge {A_*(s)\sqrt{s}\over2}
       \int_s^{2s}r^{-2}\,dr={A_*(s)\over4\sqrt{s}},
   \]

   so (2.8) has the correct endpoint and factor.  Combining this with
   `H(s)<=p`, (2.10), `J(s)=sqrt(s)H(s)`, and
   `H(s)>=c_0/sqrt(K)` gives

   \[
   {A_*(s)\over J(s)}
   \le {4\sqrt Kp\over c_0}(\beta+4u_2)
   \le {4(\beta+4u_2)\over c_0(1-\beta)},
   \]

   exactly (2.12).  There is no missing power of `alpha` in this conversion.

4. Under `beta<=sqrt(alpha)<=1/2`,
   `u_2<=2sqrt(2alpha)/c_G`.  Thus the value of `C_*` in (2.13) gives (2.14)
   exactly.  Numerically `C_*` is about `1218`, so the claim `C_*<1220` is
   correct.  Equation (2.15) follows exactly from
   `J(s)=sqrt(s)H(s)>=c_0sqrt(alpha)`.

## 2. Central-good Fisher seed

1. Equations (3.3)--(3.4) are consistent because
   `varphi(M_alpha)=I_0 exp(-M_alpha^2/2)=c_0sqrt(alpha)/8`.
   The tail estimate (3.6) has no missing factor two: on the union of the two
   tails, `I(g)<=I(delta_alpha)`, and that union has probability at most one.

2. For `eta_s>0`, (3.7) is the ordinary Markov estimate in the measure
   `q_sI(g)dy`.  If `eta_s=0`, then `tau=0` and the displayed quotient
   `A_*/tau` is `0/0`; this case must be separated.  Since then `1-e=0`
   `q_sI(g)dy`-almost everywhere, the desired conclusion still holds without
   Markov.  This is a convention gap, not a loss in the estimate.

3. The removal argument actually gives `5J(s)/8` in (3.8), so the stated
   `J(s)/2` is safe.  From `eI(g)<=I_0`, (3.9) follows with the exact constant
   `c_0/(2I_0)=1/8`.

4. The identity

   \[
   \operatorname{tr}R_{\mathcal G}
   =\int_{\mathcal G}q_s{eI(g)^2\over g(1-g)}
   \]

   and `I(g)>=c_Gg(1-g)` give (3.13).  Its last constant is exactly
   `c_Gc_0/2=1/(8pi)`.  Conditional covariance Cauchy--Schwarz gives
   `vv^T/[g(1-g)]<=Cov(X|Y)`, and total covariance gives
   `R_G<=R<=I/s`; hence (3.15) is correct.

5. **Scale error in (3.16).**  From (3.2),

   \[
   e={|v|^2\over sI(g)^2},\qquad
   \sqrt e={|v|\over\sqrt{s}\,I(g)}.
   \]

   The general angular theorem at strong-convexity parameter `t=1/s` uses

   \[
   \epsilon_y=1-{\sqrt t|v|\over I(g)}
   =1-{|v|\over\sqrt{s}\,I(g)}=1-\sqrt e.
   \]

   Thus the numerator in (3.16) must be `|v|`, divided by `sqrt(s)I(g)`;
   the displayed `sqrt(s)|v|/I(g)` is wrong by a factor `s`.  The subsequent
   use in (4.5) is consistent with the corrected formula, so no later power
   changes.

6. At the allowed endpoint `tau=1`, `G` may contain points with `v=0`, where
   `u` is undefined.  Either choose `alpha<alpha_*` so `tau<1`, or replace
   `G` by `G intersect {v!=0}`.  The latter loses no `q_sI(g)e` mass and hence
   preserves (3.8)--(3.15).

## 3. Derivative and wedge audit

1. Differentiating `v=E[fX]-gm` in the natural parameter gives

   \[
   \partial_{c_j}v_i
   =E[(f-g)(X_i-m_i)(X_j-m_j)]=D_{ij},
   \]

   so (4.2) is exact.  Since `c=y/s`,
   `nabla_yu=P_uD/(s|v|)`.  Multiplication by the amplitude
   `|v|/[sqrt(s)sqrt(g(1-g))]` proves (4.4) with the displayed
   `s^{-3/2}` scaling.

2. The posterior is `1/s`-strongly log-concave.  The general angular theorem
   therefore gives `||P_uD||_HS^2<=C_delta s^2 Omega_delta(epsilon)`.  After
   (4.4), centrality absorbs `1/[g(1-g)]`, leaving exactly `C_delta/s` in
   (4.6).

3. From `R<=I/s`, condition `s tr R>4` implies
   `||R||op/trR<1/4`.  The rank-sensitive wedge inequality then yields
   `E||nabla F||_HS^2>=trR/[2(K+s)]`.  With (3.13), this is exactly (4.8),
   and `K>32pi/alpha^(3/2)` is exactly the threshold ensuring (4.7).

4. Dividing the right side of (4.6) by the right side of (4.8) gives, up to
   the absorbed factor `16pi`, (4.9).  Since
   `Omega(O(sqrt(alpha)))=O(alpha^(1/12))` up to a logarithm, the exponent
   `-17/12` is correct.  The decomposition (4.10) is an exact orthogonal
   decomposition in the output index.

## 4. Laplace counterexample in Section 5.1

The example is correct.  For the symmetric exponential law,
`I_lambda(a)=min(a,1-a)` and `psi_lambda=1`.  Gaussian monotone likelihood
ratio makes `g_s` strictly increasing; convolution with the Gaussian kernel
then makes `T_s1_S` strictly increasing.  Every nontrivial level is a
half-line and saturates the one-dimensional profile, proving (5.7).

The Fisher feature is smooth.  If its derivative vanished `q_s`-almost
everywhere, positivity of `q_s` would make it constant everywhere; the
bounded nonconstant function `2arcsin sqrt(g_s)` would then be affine, an
impossibility.  Hence (5.9) and the refutation of (5.10) follow.  If the
ambient setup is required to be isotropic, rescale this variance-two law by
`1/sqrt(2)`; all three assertions are invariant under the corresponding
space/time scaling.

## 5. Corrections needed in the Section 6 implication

1. The word `modulus` should explicitly mean nondecreasing.  Otherwise
   `tau<=8C_*sqrt(alpha)` does not imply
   `Xi(tau)<=Xi(8C_*sqrt(alpha))`.  Equivalently, replace the second condition
   in (6.3) by
   `sup_{0<=r<=8C_*sqrt(alpha)}Xi(r)<=c_sp/4`.

2. Threshold (6.5) controls `N_sp/N`, but it need not make `N>2` if `N_sp`
   is zero or very small.  The affine contradiction separately requires

   \[
   {\alpha^{3/2}K\over8pi}>{1\over1-Xi(tau)}.
   \]

   Add this universal threshold to the final maximum.  Alternatively state
   without loss that `N_sp>=c_sp` (enlarging `N_sp` only weakens (6.2)); then
   (6.5) gives `N>4` and the written contradiction is valid.

3. `C_rad` is used as a universal constant but is not included among the
   constants asserted by the proposed lemma.  It must be added to its first
   sentence.  The parameter `N>0` should also be explicitly quantified in
   (6.1).

4. The lemma assumes merely `tr R_G>=r` for an arbitrary `r>0`, while its
   alternatives and constants contain no dependence on `r`.  The application
   only needs the fixed lower bound `r=sqrt(alpha)/(8pi)`.  As stated, however,
   the lemma also covers an arbitrarily small high-rank packet selected from
   an otherwise affine near-minimizer; neither the relative affine branch nor
   a fixed positive coarea charge is stable under such a vanishing packet.
   A credible statement must either impose `tr R_G>=r_sp` for a fixed
   universal `r_sp>0`, allow `c_sp`/`Xi` to depend quantitatively on `r`, or
   require `G` to carry a fixed fraction of the total Fisher flux.  This does
   not invalidate the algebraic implication for (6.4), but it makes the
   proposed theorem materially overstrong as currently quantified.

After repairs 1--3, and assuming a phase-charge theorem in a form that covers
the fixed trace threshold in (6.4), the remaining sufficiency calculation is
correct: (6.5) makes `N_sp/N<=c_sp/4`, the phase charge is at least
`c_sp/2`, and (5.3), (6.3) make it at most `c_sp/4`.

## 6. Disposition of audit findings

The audited source was subsequently patched as follows: (3.16) now has
`|v|/[sqrt(s)I(g)]`; the `eta_s=0` and `tau=1` cases are separated; one
fully frozen choice of `alpha`, `tau`, and `delta` is displayed in
(3.17)--(3.20); `Xi` is nondecreasing; `C_rad` and `N` are quantified; the
large-`K` threshold separately forces `N>2`; and the proposed lemma assumes
the fixed trace lower bound `sqrt(alpha)/(8pi)` used by the application.
With those changes, no correction remains in the proved fixed-scale theorem.
The phase-charge lemma itself remains explicitly labeled unproved.
