# Checkpoint 07: physical phase flux, finite splices, and global tube amplification

## 1. Status

No complete KLS proof is asserted.  The fixed-scale heat calculation now
produces a genuine, globally analytic physical boundary selector with a
fully audited numerical rank.  The scalar and selector-regularity escapes
have been removed.  The remaining statement is geometric: ranked physical
phase flux must be converted into a finite same-mass perimeter saving, or
classified as an affine/radial/product branch with bounded Poincare
constant.

The strongest new structural result is a finite-distance normal-tube
formula for exact regular isoperimetric regions.  It amplifies the local
`O(delta p^3)` curvature budget over distance `Theta(1/p)` to an exact
`O(delta p)` loss, charging support contact, focal time, and cut-locus
collision without a reach or Taylor remainder.  Its application to the
ranked heat levels still requires a quantitative near-minimizer transfer.

## 2. The physical heat selector

Let `K=C_P(mu)`, choose the fixed scale

\[
 \alpha=10^{-10},\qquad s=\alpha K,
\]

and use a balanced near-Cheeger set with boundary `p`.  The audited fixed
scale gives a good posterior flux `r_G`, the total flux `R`, and

\[
 m=\nabla F_0=P_sW,\qquad
 \Delta=E(R-|m|)<6.02\,10^{-5}p.                    \tag{2.1}
\]

The analytic selector

\[
 \omega_{an}={r_G\over R},\qquad
 M_{an}=\int\omega_{an}|\nabla F_0|\theta\theta^Td\mu
\]

is a literal coarea submeasure of the physical level boundaries and obeys

\[
 \boxed{\operatorname {tr}M_{an}>.005109p,qquad
 {\operatorname {tr}M_{an}\over\|M_{an}\|_{op}}>19.54.} \tag{2.2}
\]

Restricting to `omega_an>=1/2000` leaves trace `.004609p` and effective
rank `17.62`; imposing also `|m|/R>=1/2` leaves trace `.004489p` and rank
`17.16`.  Alternatively the floored analytic selector

\[
 \bar\omega=10^{-5}+(1-10^{-5})\omega_{an}
\]

is globally between `10^-5` and one and retains trace above `.0051p` and
effective rank above `18.8`.

The special Gaussian channel gives the exact dimension-free information
bounds

\[
 \int R\left\|\nabla{m\over R}\right\|_{HS}^2d\mu
 \le {8B_0\over s}\left(1+\log{M_s\over B_0}\right),       \tag{2.3}
\]

\[
 \int R{|\nabla\omega_{an}|^2\over
       \omega_{an}(1-\omega_{an})}d\mu
 \le {2B_0\over s}\left(\log{M_s\over B_0}+\log2\right). \tag{2.4}
\]

All logarithms are fixed functions of `alpha`, not of the dimension.  If
the cross-level normal variance is at least `8/17`, it transfers to the
smooth good-direction field as

\[
 \operatorname {Var}\big(E[v_Gv_G^T\mid F_0]\big)>.143.   \tag{2.5}
\]

An ambient capacity estimate at scale `K=s/alpha` would turn (2.5) into
only `alpha p=10^-10p` of scale-free energy, six orders below the required
`10^-4p`.  The selected-measure capacity must be `O(s)` or must come from a
finite physical competitor.

## 3. Central levels and the exact phase decomposition

For the coarea matrices `M_r`, put

\[
 a(r)=\operatorname {tr}M_r,qquad
 Q_r={M_r\over a(r)},\qquad Q={M\over\operatorname {tr}M}.
\]

Then

\[
 \boxed{
 1-\operatorname {tr}Q^2
 =E_\pi(1-\operatorname {tr}Q_r^2)
  +E_\pi\|Q_r-Q\|_{HS}^2.}                            \tag{3.1}
\]

Effective rank above seventeen forces either same-level dispersion or
between-level dispersion at least `8/17`.

Let

\[
 C=\{r:10^{-4}\le\mu(F_0>r)\le1-10^{-4}\}.
\]

The floored selector and the coarea deficit show

\[
 \operatorname {tr}M_C>.0048495p,qquad
 \operatorname {rank}_{eff}M_C>17.86,                 \tag{3.2}
\]

and, pointwise on `C`,

\[
 a(r)>1.99998\,10^{-9}p.                              \tag{3.3}
\]

For every measurable `E subset C`, the normalized level law satisfies

\[
 \pi_C(E)\le206.21|E|+.01242.                         \tag{3.4}
\]

Thus extreme levels, vanishing selected trace, and selector zero sets
cannot carry the ranked rotation.

In selected-mass coordinate, the between-level branch has the exact phase
energy

\[
 \operatorname {tr}M_C
 \int {\|Q_r'\|_{HS}^2\over a(r)}dr
 \ge {8\pi^2\over17},                                \tag{3.5}
\]

in the relaxed sense.  Charging this energy is the unresolved step.

## 4. Exact heat kinematics and the three residual sources

Write `g=|nabla F_0|`, `theta=nabla F_0/g`, `P=I-theta theta^T`, and
`H=nabla^2F_0`.  Gaussian score differentiation and the exact alignment
identity give

\[
 \boxed{\|PH\|_{op}^2\le {2I_0\over s^{3/2}}(R-g).}     \tag{4.1}
\]

Consequently the coarea-weighted physical-arclength rotation obeys

\[
 \int\omega|PH\theta|d\mu
 \le\left({2I_0\Delta\over s^{3/2}}\right)^{1/2}.       \tag{4.2}
\]

After division by the central selected flux, this tends to zero like
`alpha^-3/4 K^-1/2`.  The level-synchronized rotation instead contains the
uncontrolled lapse:

\[
 \mathcal R_{lev}=\int\omega{|PH\theta|\over g}d\mu.    \tag{4.3}
\]

On a regular band, with `Z=theta theta^T`, surface measure
`d sigma_r=omega rho dH^(k-1)`, and flow `v=theta/g`, the exact transport
equation is

\[
 \boxed{
 Q_r'={1\over a(r)}\int(v\cdot\nabla Z)d\sigma_r
 +{1\over a(r)}\int c(Z-Q_r)d\sigma_r,}              \tag{4.4}
\]

where

\[
 c={H_\Sigma-\theta\cdot\nabla V\over g}
   +{\theta\cdot\nabla\log\omega\over g}.             \tag{4.5}
\]

The first term costs `sqrt(2) R_lev`.  The two smooth source actions are at
most

\[
 \sqrt2\int\omega|H_\Sigma-\theta\cdot\nabla V|d\mu,
 \qquad \sqrt2\int|\nabla\omega|d\mu.                 \tag{4.6}
\]

The selector information controls only the flux-weighted version of the
second integral.  Critical values and hard-support contact add a matrix
Radon residual recording sheet birth, merger, and death.  The remaining
cross-level mechanisms are therefore exactly: lapse, orientation-dependent
smooth expansion/selector reweighting, and focal/contact residual.

An exact heat-generated log-affine facet model realizes all current
Gaussian equalities on the flat cores.  Its phase weights are

\[
 \pi_i(z)\propto A_i e^{-c_i^2+\sqrt2c_i z},
 \qquad c_i=\sqrt s\,\kappa_i,                         \tag{4.7}
\]

and can have vanishing same-level dispersion with between-level dispersion
`17/18`.  All charge lies at ridges, ends, or global phase transitions.
Thus Gaussian regularity alone cannot close (3.5).

## 5. Finite physical splices

For the smoothed level family, the supremum of all exactly mass-preserving
finite local splice savings satisfies

\[
 \boxed{\mathfrak G_{fin}(F_0)le6.02\,10^{-5}p.}       \tag{5.1}
\]

A planar wedge bevel with normal angle `gamma` saves, to first order, the
ridge weight times `1-cos(gamma/2)`.  For facet weights `a_i`, normals
`n_i`, admissible ridge conductances `w_ij`, and graph spectral gap
`lambda`, the simultaneous bevel saves at least

\[
 \boxed{{\lambda A\over8}\left(1-{1\over R}\right),
 \quad A=\sum_i a_i,
 \quad R={A\over\|\sum_i a_in_in_i^T\|_{op}}.}         \tag{5.2}
\]

At `A>=.004p`, `R>=17`, and `lambda>=1/4`, this exceeds
`1.17*10^-4p` and contradicts (5.1).  The numerical-expander ridge branch
is therefore closed.

An explicit simultaneous rounded box for independent one-sided exponentials
also gives a fixed improvement below the max-box perimeter.  This removes
that model as an extremizer; the coordinate halfspace remains the bounded
one-dimensional product branch.

Two exact no-go results delimit the splice theorem.

1. Moving disjoint flat log-affine patches with a common normal slope
   `kappa` cannot save perimeter.  For a graph height `h`,

   \[
   \Delta P+\kappa\Delta\mu
   =\int we^{-\kappa h}(\sqrt{1+|\nabla h|^2}-1)\ge0.  \tag{5.3}
   \]

2. Two-patch first and second variations contain no angle between spatially
   disjoint normals.  Global rank may be entirely cross-level, with every
   individual `Q_r` rank one.

Different log-affine slopes do admit a first-order volume-exchange saving.
In the ideal softmax model a between-level variance forces slope dispersion
provided the level law has controlled capacity, and an `O(1/s)` cutoff cost
turns

\[
 s\sum_i a_i(\kappa_i-\bar\kappa)^2                  \tag{5.4}
\]

into physical saving.  A multimodal level law can evade this capacity
estimate.  Once low ridge incidence has been converted into genuine
perimeter-additive components, however, an exact component charge is
available:

\[
 P(A)-\psi\mu(A)
 =\sum_i a_i\left(1-{\psi\over\kappa_i}\right).       \tag{5.5}
\]

Thus significant separated log-affine components of a near-Cheeger set
must have the common slope `kappa_i=(1+o(1))psi`; then
`sqrt(s)kappa_i=O(sqrt(alpha))` and (4.7) cannot produce fixed cross-level
rotation.  The missing physical step is low ridge/end conductance to actual
component cutting with the same quantitative error.

In the exact separated-component model the constants already close.  Slopes
`kappa_i>=2psi` carry at most

\[
 {2(6.02\,10^{-5})\over.004489}<.02683              \tag{5.6}
\]

of the selected flux.  All remaining dimensionless slopes differ by
`O(sqrt(alpha))`.  Centrality of `mu(A_r)` does **not** bound the heat
profile coordinate, so the first bounded-interval likelihood-ratio argument
was rejected in audit.  The repaired unbounded-coordinate lemma uses the
second moment of the shifted-Gaussian phase mixture and the selector floor:

\[
 \sqrt{\operatorname {Var}Q}
 \le {\Delta_c\over\sqrt{2\omega_0}}
       \sqrt{1+2C_c^2}+2\sqrt q.                     \tag{5.7}
\]

It proves that cross-level variance `8/17` forces cut-plus-tail error larger
than `1.4*10^-4p`.  No Poincare inequality for the multimodal level law or
bounded profile interval is used.

The first attempt to identify small exterior ridge conductance with a cheap
component cut was rejected.  If exterior phase traces `U,V` on
`partial E` are separated by an interior finite-perimeter cut, the exact BV
identity is

\[
 P(E_1)+P(E_2)=P(E)+2\operatorname {Fill}_E(\partial U). \tag{5.8}
\]

A codimension-two ridge tube does not control this filling.  A long
rectangle and a ball give explicit arbitrarily large ridge-to-fill ratios.
The low-ridge branch therefore has a necessary third alternative:

\[
 \text{cheap interior fill/component synchronization}
 \quad\text{or}\quad
 \text{large filling/coherent or concurrent organization}.    \tag{5.9}
\]

The second branch is being coupled to the long-ray covariance and translated
thin-shell mechanisms; omitting it would make the component argument false.

## 6. Profile linearity and exact equality rigidity

If a union of separated components has volumes `v_i`, then

\[
 P(\cup_iA_i)-I(v)
 \ge\sum_i v_i\{q(v_i)-q(v)\},
 \qquad q(t)=I(t)/t.                                  \tag{6.1}
\]

Exact separated equality forces `I(t)=q(v)t` on `[0,v]`, not merely at the
component masses.  Quantitatively, if
`q(v_i)<=(1+delta)q(v)` and `v_i<=(1-eta)v`, then

\[
 q(t)\le q(v)\left(1+{\delta\over\theta\eta}\right)
 \quad(\theta v_i\le t\le v_i).                       \tag{6.2}
\]

In one dimension

\[
 -{d\log q(v)\over d\log v}
 =1-{I'(v)\over q(v)}.                                \tag{6.3}
\]

Thus exact linearity is an exponential tail, and near linearity gives an
explicit conditional-tail approximation to an exponential law.

For smooth full support, exact profile linearity makes every regular
minimizer a hyperplane; distinct level hyperplanes are parallel.  The
problem descends to a variance-one one-dimensional log-concave marginal,
whose Cheeger constant is universal.  The same exact classification now
holds for regular hard support.  The free-boundary second variation includes
the nonnegative contact term

\[
 \int_{\partial\Sigma}II_{\partial\Omega}(N,N),        \tag{6.4}
\]

and two nonparallel planar free-boundary sections of a `C^1` convex body
intersect in its interior.  They cannot be hidden in the codimension-seven
singular set.  Hence the exact zero-curvature branch again reduces to one
affine marginal.

The quantitative local curvature budget alone does not close: a thin neck
in hypersurface dimension at least four can rotate a fixed normal angle
while `int|II|^2` tends to zero.

## 7. Global killed-tube amplification

Let `E` be an exact regular isoperimetric region of volume `v_0`, perimeter
`P_0=I(v_0)`, and constant weighted mean curvature `lambda`.  Along an
exterior normal ray, stop at first support contact, focal time, or loss of a
unique nearest point.  Before killing, the exact weighted Jacobian is

\[
 j_x(t)=e^{\lambda t-D_x(t)},                           \tag{7.1}
\]

\[
 D_x(t)=\int_0^t(t-u)\left[
 \operatorname {tr}\{S_x^2(I+uS_x)^{-2}\}
 +\nabla^2V(x+uN_x)[N_x,N_x]\right]du\ge0.             \tag{7.2}
\]

The killed normalized flux

\[
 R(t)=\int_\Sigma1_{\{t<\tau(x)\}}e^{-D_x(t)}d\sigma_\mu(x)
\]

is nonincreasing.  Suppose `0<a<v_0<b<=1/2`,

\[
 c=I(b)/b,\quad \beta=b/v_0,\quad
 I(a)/a\le(1+\delta)c,
\]

and set `s_0=min(v_0-a,b-v_0)`, `eta=delta b/s_0<1`.  If `T_b` is the first
outer-tube time reaching volume `b`, then

\[
 {1\over(1+\eta)c}\log\left(1+{1+\eta\over1+\delta}
 (\beta-1)\right)
 \le T_b\le{\log\beta\over c},                       \tag{7.3}
\]

and

\[
 \boxed{P_0-R(T_b-)\le
 \left(1-{\beta^{-\eta}\over1+\delta}\right)P_0
 \le(\delta+\eta\log\beta)P_0.}                     \tag{7.4}
\]

This is the desired `p^3 times (1/p)^2=p` amplification with exact singular
accounting.

If a packet `G` survives to `T` with `D_x(T)<=h`, then

\[
 \boxed{
 \operatorname {Cov}(\mu)\succeq
 {e^{-|\lambda|T-h}T^3\over12}
 \int_GN\otimes N\,d\sigma_\mu.}                    \tag{7.5}
\]

This closes a coherent-normal packet in isotropic position.  It does not
close arbitrarily high normal rank.  Independent exponentials with a union
of rare coordinate tails have exact normal matrix `(P/m)I`, ridge/perimeter
ratio asymptotic to `alpha/2`, and killed-flux loss only `Theta(alpha)P`.
The model calibrates the required collision/conductance-versus-product
inverse theorem.

Two gaps precede any use of (7.4)--(7.5) in the heat proof:

1. the ranked physical packet lies on integrated near-minimizer levels of
   `F_0`, not on exact isoperimetric minimizers;
2. in the high-rank survivor, small collision loss must imply either a
   bounded-reuse ridge saving or a quantitative almost-product splitting.

Neither transfer may be hidden under regularity or compactness.

## 8. Locked-line and random-fiber audits

For hard mass-preserving stochastic localization, an exactly fixed protected
line forces the original set to have conditional mass one half on almost
every parallel line.  Gaussian convolution injectivity supplies this
statement.  If `V_u` is the mean conditional variance along the line, then

\[
 P_\mu(S)\ge{1\over\sqrt{48V_u}}.                     \tag{8.1}
\]

For a deterministic direction in isotropic position, `V_u<=1`, so this
branch is closed.  Approximate locking gives an explicit finite-time bound.
The unresolved quantity for a history-selected line is

\[
 E{p_\tau^{3/2}\over\sqrt{K_\tau}},                   \tag{8.2}
\]

which isotropy alone does not control.  A bounded-speed inertial tracker has
an exact entropy/kinetic identity, but calibrated entropy bursts can still
move the active line without a universal good-event probability.

A separate Haar-random fiber descent is blocked sharply: on the cube,
simplex, and shifted exponential product, conditional fiber variance is
`Theta(1/n)` while every linear mode receives only `Theta(1/n^2)` of its
variance as fiber energy.

## 9. Active and blocked mechanisms

### Active

1. Prove a low-incidence phase cut: small ridge/end/contact conductance must
   produce genuine perimeter-additive components with controlled error.
2. Use (5.5) to synchronize the slopes of all significant separated
   log-affine components, eliminating fixed cross-level softmax rotation.
3. Prove the high-rank killed-tube inverse: small collision loss implies an
   almost-product decomposition whose one-dimensional factors have bounded
   variance-normalized Cheeger constants.
4. Transfer the ranked near-minimizer heat levels to a quasiminimal or exact
   tube setting without losing the normal matrix.

### Blocked

- Gaussian selector information without a selected-measure capacity theorem;
- any direct bound of level rotation by physical-arclength rotation, because
  of the exact factor `1/|nabla F_0|`;
- local `L^2` curvature-to-angle conversion in high hypersurface dimension;
- translations of common-slope flat components;
- generic random fibers, scalar Bernstein propagation, and an ambient
  Poincare estimate at scale `K`;
- calling a low-ridge high-rank survivor “product” without proving the
  quantitative splitting and covariance transfer.

## 10. Audit ledger

- All new reports render through Pandoc.
- The analytic-selector constants use the conservative certified lower
  bound `b_0>.00517`.
- The same-level/cross-level matrix identity and every numerical rank loss
  were recalculated after central restriction.
- The global tube theorem is under an independent clean-room audit.
- The tube theorem has not been applied to the heat levels; the exact-versus-
  near-minimizer distinction is explicit.
- The product-exponential max-tail family is used only as a stress test, not
  as an isoperimetric minimizer.
- No dimension-free Poincare, selected boundary capacity, or KLS-equivalent
  covariance-process estimate is assumed.
