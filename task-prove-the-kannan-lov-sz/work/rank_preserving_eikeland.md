# Rank-preserving proximal minimization: what is true and what is false

## 0. Outcome

Let `A` be a finite-perimeter set of volume `v` and write

\[
             \delta=P_\mu(A)-I_\mu(v)\ge0.                 \tag{0.1}
\]

There is a completely dimension-free proximal selection.  For every
`Lambda>0`, a minimizer `B` of

\[
 P_\mu(E)+\Lambda\mu(E\mathbin\triangle A),
             \qquad \mu(E)=v,                              \tag{0.2}
\]

satisfies

\[
 \boxed{\quad P_\mu(B)\le P_\mu(A),\qquad
 \mu(A\mathbin\triangle B)\le {\delta\over\Lambda},\qquad
 P_\mu(B)\le P_\mu(C)+\Lambda\mu(B\mathbin\triangle C)\quad} \tag{0.3}
\]

for every same-volume `C`.  On every regular part of its boundary this
implies

\[
                   |H_\mu-\lambda|\le\Lambda              \tag{0.4}
\]

for a constant `lambda`.  Thus the scalar Ekeland step itself loses no
dimension and gives exactly the error which the killed-tube formula can
tolerate: over a tube of length `T`, the constant-CMC formulas acquire only
the factors `exp(+-Lambda T)`.

If, in addition, the proximal problem admits a bounded calibration `z` with

\[
 |z|\le1,\quad z\cdot n_B=1,quad
 \operatorname {div}_\mu z=\lambda+e,\quad |e|\le\Lambda, \tag{0.5}
\]

then the old boundary is sharply aligned with the calibration:

\[
 \boxed{\int_{\partial^*A}(1-z\cdot n_A)dP_\mu\le2\delta.} \tag{0.6}
\]

The potentially large multiplier `lambda` cancels because `A` and `B` have
the same volume.  This is stronger than a bound containing
`|lambda| mu(A triangle B)`.

However, (0.3)--(0.6) do **not** transfer the quadratic normal matrix from
`partial A` to `partial B`.  Divergence controls scalar flux; transporting
`z tensor z` requires a curl or derivative estimate for `z`, and neither
Ekeland's principle nor (0.5) supplies one.  A local divergence-free rotating
field below makes this obstruction exact.

There is a variational way to preserve the matrix: add its nuclear-norm
distance to the objective.  For `0<kappa<1/3`, this gives a lower
semicontinuous uniformly elliptic problem and the exact estimate

\[
 \boxed{\|M_\omega(B)-M_\omega(A)\|_*\le\delta/\kappa.}   \tag{0.7}
\]

Its Euler equation is anisotropic and contains `kappa nabla omega`.  The
available Gaussian-selector regularity does not make that term perturbative
on a tube of length `asymp1/p`; the fixed-scale numerical audit is given in
Section 6.  Consequently (0.7) is not a hidden completion of the transfer.

Finally, an explicit isotropic one-dimensional asymmetric Laplace measure
shows that replacement by an exact isoperimetric minimizer cannot preserve
even the trace of a positive analytic selected packet.  At the numerical
values of the fixed-scale argument, the starting set has deficit below
`6.02*10^-5 p`, selected flux exactly `.005p`, while the unique exact
minimizer retains less than `2*10^-5` of that packet.  Thus any valid transfer
theorem must genuinely use the quantitative-stationarity branch; exact
minimizer replacement and ordinary compactness are false.

## 1. The proximal BV lemma

Work on the affine support `E` of `mu`, with its Euclidean metric.  Relative
weighted `BV` perimeter is denoted by `P_mu`.  No smoothness is needed in
this section.

**Lemma 1.1 (dimension-free proximal selection).**  Let `mu` be a
log-concave probability on `E`, not a point mass.  Let `A` have finite
relative weighted perimeter and `0<v=mu(A)<1`.  For every `Lambda>0`, problem
(0.2) has a minimizer.  Every minimizer satisfies (0.3).

**Proof.**  A minimizing sequence has uniformly bounded weighted perimeter.
On every compact subset of the relative interior of the support, the density
is bounded below by a positive number, so ordinary local `BV` compactness
applies.  A diagonal subsequence converges locally in `L^1`.  Since `mu` is
a probability, choose compact sets whose complements have arbitrarily small
`mu`-mass; local convergence is therefore global `L^1(mu)` convergence.
Volume and the fidelity term pass to the limit, and weighted relative
perimeter is lower semicontinuous.  This proves existence, including on an
unbounded support.  It does not assume attainment of `I_mu(v)`.

Comparison with `A` gives

\[
       P_\mu(B)+\Lambda\mu(A\triangle B)\le P_\mu(A).
\]

Since `P_mu(B)>=I_mu(v)=P_mu(A)-delta`, the first two conclusions in (0.3)
follow.  For a same-volume `C`, minimality and the triangle inequality give

\[
\begin{aligned}
 P_\mu(B)+\Lambda\mu(B\triangle A)
 &\le P_\mu(C)+\Lambda\mu(C\triangle A)\\
 &\le P_\mu(C)+\Lambda\mu(C\triangle B)
              +\Lambda\mu(B\triangle A),
\end{aligned}
\]

which is the last assertion.  QED.

The same conclusions follow from Ekeland's variational principle without
attainment of (0.2), but the direct minimizer is stronger and costs nothing.

## 2. Quantitative stationarity and the perturbed tube identity

Assume temporarily that the density is `C^2` and positive on a smooth
support and that a relatively open piece `Sigma` of `partial^*B` is `C^2`.
For every smooth volume-preserving normal speed `u`, (0.3), applied to both
signs of the flow, yields

\[
 \left|\int_\Sigma H_\mu u\,d\sigma_\mu\right|
       \le\Lambda\int_\Sigma|u|\,d\sigma_\mu,
 \qquad \int_\Sigma u\,d\sigma_\mu=0.              \tag{2.1}
\]

The quotient norm identity

\[
 \inf_{c\in\mathbb R}\|H_\mu-c\|_{L^\infty}
 =\sup_{\int u=0}{|\int H_\mu u|\over\int|u|}       \tag{2.2}
\]

then proves (0.4).  Equivalently, the essential oscillation of `H_mu` is at
most `2Lambda`, and `lambda` can be chosen as the midpoint of its essential
range.  On a smooth convex hard support the regular free boundary meets the
support orthogonally: after blow-up the fidelity term is lower order, so the
usual relative-perimeter halfspace condition is unchanged.

Write `H_mu(x)=lambda+e_x`, `|e_x|<=Lambda`.  Along every surviving exterior
normal ray the exact Jacobian is

\[
 j_x(t)=\exp\{(\lambda+e_x)t-D_x(t)\},\qquad D_x(t)\ge0. \tag{2.3}
\]

With `tau(x)` the first focal/contact/cut time and

\[
 R(t)=\int_\Sigma1_{\{t<\tau(x)\}}e^{-D_x(t)}d\sigma_\mu(x), \tag{2.4}
\]

one has, at every regular tube time,

\[
 e^{(\lambda-\Lambda)t}R(t)
 \le P_\mu(B_t)\le e^{(\lambda+\Lambda)t}R(t).       \tag{2.5}
\]

Thus a constant-CMC killed-tube calculation changes by at most
`e^{Lambda T}` up to time `T`.  In particular, choosing
`Lambda=eta q` when the relevant profile slope is `q` and `T<=C/q` costs
only `e^{C eta}`.  Formula (2.5) is literal, not a Taylor expansion.

There are two regularity qualifications which cannot be omitted.  A bounded
mean-curvature quasiminimizer is smooth only off its singular set, which can
be nonempty in dimension at least eight.  Also `C^{1,alpha}` regularity by
itself does not give positive reach.  Therefore use of (2.5) in arbitrary
dimension still requires the rectifiable killed-tube treatment of the
singular projection; (0.3) does not make that issue disappear.

## 3. What a calibration would give

Assume (0.5), with normal traces understood in the theory of
divergence-measure fields.  Weighted Gauss--Green and `mu(A)=mu(B)` give

\[
\begin{aligned}
 \int_{\partial^*A}z\cdot n_A\,dP_\mu-P_\mu(B)
 &=\int(1_A-1_B)\operatorname {div}_\mu z\,d\mu\\
 &=\int(1_A-1_B)e\,d\mu.                              \tag{3.1}
\end{aligned}
\]

Consequently

\[
\begin{aligned}
 0\le\int_{\partial^*A}(1-z\cdot n_A)dP_\mu
 &\le P_\mu(A)-P_\mu(B)+\Lambda\mu(A\triangle B)\\
 &\le2\delta,                                         \tag{3.2}
\end{aligned}
\]

which proves (0.6).  In particular,

\[
 \int_{\partial^*A}|n_A-z|^2dP_\mu\le4\delta.        \tag{3.3}
\]

For any selector `0<=omega<=1`, putting
`P_omega(A)=int omega dP_mu`, one gets the explicit same-surface projector
bound

\[
 \int_{\partial^*A}\omega
   \|n_A n_A^T-zz^T\|_{op}dP_\mu
 \le4\sqrt{P_\omega(A)\delta}.                       \tag{3.4}
\]

This is the full conclusion of scalar calibration.  It compares `n_A` with
`z` **on `partial A`**, not with `n_B` on `partial B`.

Existence of (0.5) is stronger than (0.3).  The fixed-volume set problem is
not a convex problem: its value as a function of volume need not have a
supporting affine functional.  A max-flow/min-cut calibration is available
for the corresponding unconstrained convex `TV` problem, but it cannot be
inserted here without proving exact selection of the prescribed volume.
Thus (0.5) is a conditional strengthening, not an automatic consequence of
Lemma 1.1.

## 4. Scalar divergence cannot transport the normal matrix

To compare the tensor on `partial A` with that on `partial B`, Gauss--Green
would have to be applied to

\[
                    (z^THz)z.                         \tag{4.1}
\]

Its weighted divergence contains

\[
 2\langle Hz,(\nabla z)z\rangle+(z^THz)\operatorname {div}_\mu z. \tag{4.2}
\]

The first term is uncontrolled by (0.5).  It is a curl/rotation term, not a
technical artifact.

Here is an exact local model.  In the plane, away from the origin, let

\[
 z(x)={1\over|x|}(-x_2,x_1)=e_\theta(x).             \tag{4.3}
\]

Then `|z|=1` and `div z=0`.  Every radial line segment is a hypersurface
whose unit normal is `+-z`.  Choose a convex ball which does not contain the
origin but intersects two radial rays making a fixed angle `vartheta>0`.
The two ray segments are disjoint regular hypersurface patches, scalar
calibration is saturated on both, and the divergence error is zero, while

\[
 \|n_1n_1^T-n_2n_2^T\|_{op}=|\sin\vartheta|.          \tag{4.4}
\]

The rays meet only at the excluded origin.  Hence neither convexity of the
ambient patch nor bounded scalar divergence prevents a fixed projector
rotation.  Any estimate transporting `zz^T` must contain a derivative/curl
quantity or a genuine geometric incidence term.

## 5. A matrix-fidelity problem which really preserves rank

Let `omega` be continuous with `0<=omega<=1`, and define

\[
 M_\omega(E)=\int_{\partial^*E}\omega(x)n_E n_E^T\,dP_\mu. \tag{5.1}
\]

For `0<kappa<1/3`, minimize at volume `v`

\[
       J_\kappa(E)=P_\mu(E)+
          \kappa\|M_\omega(E)-M_\omega(A)\|_*.       \tag{5.2}
\]

This problem has a minimizer.  The nontrivial point is lower
semicontinuity.  By nuclear/operator duality,

\[
 J_\kappa(E)=\sup_{H=H^T,\ \|H\|_{op}\le1}
 \left\{\int_{\partial^*E}
 [1+\kappa\omega(x)n^THn]dP_\mu
 -\kappa\operatorname {tr}(HM_\omega(A))\right\}.   \tag{5.3}
\]

For fixed `H`, the one-homogeneous integrand is

\[
       \Phi_H(x,\xi)=\rho(x)\left(|\xi|+
       \kappa\omega(x){\xi^TH\xi\over|\xi|}\right). \tag{5.4}
\]

If `h` is tangent to the unit vector `n`, its second derivative in the
`h` direction is

\[
 |h|^2+\kappa\omega
       \{2h^THh-(n^THn)|h|^2\}\ge(1-3\kappa)|h|^2. \tag{5.5}
\]

Thus (5.4) is convex and uniformly elliptic.  Each functional in (5.3) is
lower semicontinuous, and so is their supremum.  The direct-method argument
of Lemma 1.1 gives existence, also on an unbounded support.

Comparison with `A` and the definition of `I_mu(v)` give exactly

\[
 \kappa\|M_\omega(B)-M_\omega(A)\|_*
 \le P_\mu(A)-P_\mu(B)\le\delta,                    \tag{5.6}
\]

which is (0.7).  If

\[
 t=\operatorname {tr}M_\omega(A),\qquad
 o=\|M_\omega(A)\|_{op},\qquad e=\delta/\kappa,
\]

then

\[
 \operatorname {tr}M_\omega(B)\ge t-e,qquad
 \|M_\omega(B)\|_{op}\le o+e.                     \tag{5.7}
\]

In particular, effective rank above seventeen is retained whenever

\[
                         e<{t-17o\over18}.           \tag{5.8}
\]

On the regular part, nonsmooth variational calculus supplies a symmetric
subgradient `H`, `||H||op<=1`.  The Euler equation is the prescribed
anisotropic-mean-curvature equation for

\[
                    \phi(x,n)=1+\kappa\omega(x)n^THn. \tag{5.9}
\]

It contains both `kappa omega II` and `kappa nabla omega`.  It is not the
Euclidean CMC equation (0.4).

## 6. Fixed-scale numerical and derivative audit

For the floored analytic packet in
`fixed_scale_physical_splicing.md`, one only knows

\[
 t>.0051p,\qquad t/o>18.8,
 \qquad \delta_{co}\le6.02\,10^{-5}p.               \tag{6.1}
\]

The rank margin in (5.8), using only these audited values, is

\[
 {t-17o\over18}
 >{.0051\over18}\left(1-{17\over18.8}\right)p
 >2.71\,10^{-5}p.                                   \tag{6.2}
\]

Since uniform ellipticity requires `kappa<1/3`, (5.6) cannot reach (6.2)
with the current deficit: even at the endpoint it gives an error larger
than `1.80*10^-4p`.  Retuning the fixed heat fraction can reduce the scalar
deficit and make (5.8) true, but it does not solve the Euler-equation issue.

Indeed the Gaussian information estimate controls spatial selector
variation only at the scale

\[
                         |\nabla\omega|\sim s^{-1/2}. \tag{6.3}
\]

The relevant tube has length `T asymp1/p`, while `s=alpha K` and the
worst KLS scaling is `p sqrt K asymp1`.  Therefore the available bound on
the accumulated anisotropic forcing is of size

\[
                  \kappa T/\sqrt s\asymp
                  {\kappa\over\sqrt\alpha}.         \tag{6.4}
\]

On the other hand, the coarea deficit is `delta_co/p=O(sqrt alpha)`.
Preserving the rank margin `gamma p` through (5.6) requires
`kappa >= c sqrt(alpha)/gamma`.  Substitution into (6.4) leaves the fixed
loss `c/gamma`; here `gamma` is only of order `10^-5`--`10^-4`.  Thus the
currently proved derivative estimate gives a loss of order at least
`10^4`, independently of how small `alpha` is.  This is the same fixed-scale
power obstruction in a different variational form.

This calculation does not say that every anisotropic minimizer realizes the
upper bound (6.4).  It says that the audited selector regularity is
insufficient to justify treating (5.9) as the Euclidean CMC equation.

## 7. Exact-minimizer replacement is false

Fix `tau>0` and consider the asymmetric Laplace density on the line

\[
 f_\tau(x)=c_\tau
 \begin{cases}
 e^x,&x\le0,\\
 e^{-(1+\tau)x},&x\ge0,
 \end{cases}
 \qquad c_\tau={1+\tau\over2+\tau}.                 \tag{7.1}
\]

It is log-concave.  For
`0<v<1/(2+tau)`, the left and right half-lines of mass `v` have perimeters

\[
                         P_L(v)=v,qquad
                         P_R(v)=(1+\tau)v.           \tag{7.2}
\]

The one-dimensional log-concave isoperimetric theorem says that half-lines
minimize; hence

\[
                         I_\mu(v)=v,                 \tag{7.3}
\]

and the left half-line is the unique minimizer modulo null sets.  The median
lies on the left exponential branch and its boundary density is

\[
                              p={1\over2}.            \tag{7.4}
\]

Take

\[
 a_0=.005,qquad \tau=.012,qquad
 v={a_0p\over1+\tau},                                \tag{7.5}
\]

and let `A` be the right half-line of mass `v`.  Then

\[
 P_\mu(A)=a_0p,qquad
 P_\mu(A)-I_\mu(v)={a_0p\tau\over1+\tau}
 <5.93\,10^{-5}p<6.02\,10^{-5}p.                   \tag{7.6}
\]

Let `q_R` be its boundary point and let `q_L` be the boundary point of the
unique left minimizer.  For `eta=10^-5`, choose `sigma>0` so small that

\[
 \omega(x)=\eta+(1-\eta)
      \exp\{-|x-q_R|^2/\sigma^2\}                   \tag{7.7}
\]

satisfies `omega(q_L)<2eta`.  This selector is positive, real analytic, and
lies between `eta` and one.  The selected boundary flux of `A` is exactly
`.005p`, whereas that of the unique exact minimizer is at most

\[
             2\eta v<{2\eta\over1+\tau}P_\mu(A).    \tag{7.8}
\]

Thus exact replacement loses more than `99.998%` of the selected packet.
After translating and scaling (7.1) to variance one, all perimeters
`p,P(A),delta` acquire the same factor, and (7.5)--(7.8) are unchanged.
The counterexample is therefore isotropic.

There is a sharp quantitative-stationarity version.  Split a small total
volume between a left tail of mass `w_L` and a right tail of mass `w_R`.
The perimeter and excess are

\[
 P=w_L+(1+\tau)w_R,qquad
 P-I(w_L+w_R)=\tau w_R.                              \tag{7.9}
\]

The two boundary components have weighted mean curvatures `1` and
`1+tau`; hence the smallest possible error in (0.4) for this set is exactly
`tau/2`.  Since its selected right flux is `(1+tau)w_R`,

\[
 {\tau\over2}={1+\tau\over2}
 {P-I\over (1+\tau)w_R}.                             \tag{7.10}
\]

This shows that a quasistationarity loss proportional to
`deficit/selected_flux` is unavoidable.  It also explains why the
quasistationary alternative is the only viable one: the exact minimizer
loses the packet, while the original two-tail set is quantitatively CMC at
the sharp scale (7.10).

## 8. Hard support, singular sets, and nonattainment audit

1. **Lower-dimensional support.**  Every statement above is made on the
   affine support.  The Euclidean metric, `BV` derivative, normal, and
   divergence are intrinsic to that support.  A point mass is excluded.
2. **Unbounded support and nonattainment of `I`.**  Lemma 1.1 and the
   matrix-fidelity problem use direct compactness of their own objectives;
   neither assumes an isoperimetric minimizer.  The inequality
   `P(B)>=I(v)` uses only the definition of the infimum.
3. **Hard support.**  Perimeter is relative perimeter.  The direct method
   remains valid for an extended-valued convex potential.  On a smooth
   support, regular proximal boundaries meet it orthogonally; on a
   nonsmooth support the correct statement is a normal-cone condition.
   No classical contact angle may be asserted at a corner.
4. **Singular set.**  Proximal and anisotropic minimizers can have a
   codimension-eight singular set.  Neither (0.3) nor uniform ellipticity
   proves that its metric projection has zero tube contribution in the
   exact generality needed by the global killed-tube theorem.
5. **Calibration.**  Equation (0.6) is fully rigorous if (0.5) exists, but
   fixed-volume Ekeland minimization alone does not prove such a global
   bounded calibration.  Local first variation proves (0.4), not (0.5).

## 9. Precise conclusion

The scalar part of the desired transfer is available with universal
constants: (0.3), (0.4), and (2.5).  Exact-minimizer replacement is
disproved by (7.1)--(7.8).  Even a hypothetical scalar calibration transfers
only alignment on the old boundary, as (4.3)--(4.4) show.  Matrix fidelity
can be imposed variationally and gives the exact rank estimate (5.6)--(5.8),
but its anisotropic Euler equation is not covered by the Euclidean killed
tube theorem and the available selector derivative estimate fails the
numerical perturbative audit (6.4).

Therefore a valid completion must add one genuinely new ingredient: either
a curl/rotation estimate for the proximal calibration, an isotropic
rank-preserving proximal functional whose first variation contains no
spatial-selector derivative, or a killed-tube theorem for (5.9) with a
finite defect bound that beats (6.4).  None of these follows from scalar
near-minimality.
