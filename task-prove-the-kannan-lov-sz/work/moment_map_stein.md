# The moment-map Stein-kernel route to KLS

## 1. Verdict

The moment map gives an exact and attractive weighted Poincare certificate.  If
`mu` is centered and isotropic, there is a positive-semidefinite matrix field
`tau` such that

\[
 \mathbb E_\mu\tau=I,
 \qquad
 \mathbb E_\mu[Xf(X)]=\mathbb E_\mu[\tau(X)\nabla f(X)],
 \qquad
 \operatorname {Var}_\mu f
 \leq \mathbb E_\mu\langle\tau\nabla f,\nabla f\rangle.       \tag{1.1}
\]

These statements are exact, including after the usual approximation for a
convex body.  They do **not** currently yield a dimension-free unweighted
Poincare inequality.  There are three distinct obstructions.

1. Comparing the two Dirichlet forms for every test function is equivalent to
   the pointwise bound `tau <= C I`.  This is false already for a one-dimensional
   shifted exponential.
2. Trace concentration, determinant control from Monge--Ampere, the
   dimension-free thin-shell theorem, and the slicing theorem do not prevent a
   logarithmic top eigenvalue.  The exact product-exponential moment map has
   all of these scalar controls and nevertheless
   `E ||tau||_op = H_d ~ log d`.
3. The correlation with the gradient is real, rather than an artifact of
   taking an operator norm.  For the product exponential, the 1-Lipschitz
   function `max_i Y_i` has unweighted energy one and bounded variance, but its
   moment-map energy is exactly `H_d`.  The isotropic simplex has the same
   phenomenon.

The correct Hardy compensation in the product model acts on the variance
*before* the weighted Brascamp--Lieb upper bound is introduced: the
one-dimensional exponential has Poincare constant at most `4`, and this
tensorizes.  No bound on the large-`tau` part of the weighted Dirichlet form can
recover that argument.

There is one logically viable residual target.  It is enough to prove a tail
estimate only for low spectral modes (or first-eigenfunction approximants) of
the ordinary Dirichlet form.  In a useful normalization it would say

\[
 \mathbb E\langle (\tau-LI)_+\nabla f,\nabla f\rangle
 \leq \theta\operatorname {Var}f+B\mathbb E|\nabla f|^2,
 \qquad \theta<1,                                      \tag{1.2}
\]

for low modes `f`, with universal `L,B,theta`.  Such a statement would close
KLS, but it is not supplied by thin shell or slicing and is not proved here.
For arbitrary smooth `f`, (1.2) is false even in one dimension.

## 2. Exact moment-measure and Stein identities

Work first in the smooth, strictly convex, full-dimensional setting.  Let
`mu` be a centered log-concave probability on `R^d`.  The moment-measure
theorem supplies an essentially continuous convex function `phi`, unique up
to translation of its argument, normalized so that

\[
 d\nu(x)=e^{-\phi(x)}dx,
 \qquad \int e^{-\phi}=1,
 \qquad (\nabla\phi)_\#\nu=\mu.                       \tag{2.1}
\]

Put

\[
 A(x)=\nabla^2\phi(x),\qquad Y=\nabla\phi(X),\quad X\sim\nu.
\]

When the gradient map is invertible, define

\[
 \tau(y)=A((\nabla\phi)^{-1}(y)).                      \tag{2.2}
\]

In the general smooth case the right definition is the conditional matrix

\[
 \tau(y)=\mathbb E_\nu[A(X)\mid\nabla\phi(X)=y].       \tag{2.3}
\]

It is positive semidefinite.  Integration by parts gives, for a smooth scalar
test `f`,

\[
\begin{aligned}
 \mathbb E_\mu[Y_i f(Y)]
 &=\int \partial_i\phi(x)f(\nabla\phi(x))e^{-\phi(x)}dx\\
 &=\int \partial_i(f\circ\nabla\phi)(x)e^{-\phi(x)}dx\\
 &=\mathbb E_\mu\sum_j\tau_{ij}(Y)\partial_jf(Y).     \tag{2.4}
\end{aligned}
\]

For vector fields `F`, equivalently,

\[
 \mathbb E\langle Y,F(Y)\rangle
 =\mathbb E\operatorname {Tr}(\tau(Y)\nabla F(Y)).    \tag{2.5}
\]

Taking `f(y)=y_j` shows

\[
 \mathbb E_\mu\tau=\operatorname {Cov}(\mu).          \tag{2.6}
\]

In particular, `E tau=I` in isotropic position.

Apply Brascamp--Lieb to `g=f\circ\nabla\phi` under `nu`.  Since
`nabla g=A nabla f(nabla phi)`,

\[
\begin{aligned}
 \operatorname {Var}_\mu f
 &=\operatorname {Var}_\nu(f\circ\nabla\phi)\\
 &\leq \mathbb E_\nu
   \langle A^{-1}\nabla(f\circ\nabla\phi),
                 \nabla(f\circ\nabla\phi)\rangle\\
 &=\mathbb E_\mu\langle\tau\nabla f,\nabla f\rangle.    \tag{2.7}
\end{aligned}
\]

For nonsmooth targets, one first regularizes the convex target potential and,
if necessary, adds a small quadratic term.  The pushforward, Stein identity,
and quadratic-form inequality then pass to the limit for compactly supported
smooth tests.  In all of the explicit models below the pushforward can instead
be checked directly, so no regularity assertion is hidden in the computation.

### 2.1 Monge--Ampere and entropy

If the target density is written in normalized form `dmu(y)=e^{-V(y)}dy`,
then change of variables in (2.1) gives

\[
 e^{-\phi(x)}=e^{-V(\nabla\phi(x))}\det A(x),
 \qquad
 \log\det\tau(Y)=V(Y)-\phi(X).                         \tag{2.8}
\]

Consequently

\[
 \mathbb E_\mu\log\det\tau=h(\mu)-h(\nu).            \tag{2.9}
\]

Also, by concavity of `log det` and (2.6),

\[
 \mathbb E\log\det\tau\leq\log\det\mathbb E\tau=0. \tag{2.10}
\]

Equality forces `tau=I` almost surely in the nondegenerate smooth setting and
hence a quadratic moment potential.  The important point for the present
route is that a lower bound of order `-Cd` in (2.9), even if available, gives
no operator bound; Sections 5.3 and 5.4 give realized counterexamples.

### 2.2 Exact contractions seen by thin shell

Let `P` be an orthogonal projection and put `q(y)=y^T P y`.  Use (2.5) with
`F(y)=Py q(y)`.  Its Jacobian is

\[
 \nabla F=qP+2(Py)(Py)^T,
\]

and therefore

\[
 \mathbb E q^2
 =\mathbb E\left[q\operatorname {Tr}(P\tau)
       +2(Py)^T\tau(Py)\right].                       \tag{2.11}
\]

For `P=I`,

\[
 \mathbb E|Y|^4
 =\mathbb E\left[|Y|^2\operatorname {Tr}\tau
       +2Y^T\tau Y\right].                            \tag{2.12}
\]

The dimension-free thin-shell theorem applied to every marginal controls
`Var(q)` by a universal multiple of `rank(P)`.  Equations (2.11)--(2.12) are
the exact information this gives about the moment kernel.  It controls
particular contractions in which `tau` is coupled to `Y`; it does not control
`Tr(P tau)` unweighted, much less `||tau||_op`.  The product exponential below
satisfies all these identities while its top eigendirection is selected by a
nonlinear gradient.

The slicing theorem controls the height, and hence the entropy up to an
additive `O(d)`, of the target density.  In (2.9) the source entropy remains,
and even two-sided `O(d)` control of the entropy difference would still be
consistent with the product and simplex spectra below.

## 3. Why direct removal of the matrix weight is pointwise operator control

Here is a useful exact distinction between KLS and an overly strong proposed
intermediate statement.

**Lemma 3.1 (multiplication-form comparison).**  Let `rho` be positive on an
open set `U`, and let `T:U -> S_+^d` be locally integrable.  Then

\[
 \int\langle T\nabla f,\nabla f\rangle\,\rho
 \leq C\int|\nabla f|^2\rho
 \quad\hbox{for every } f\in C_c^\infty(U)             \tag{3.1}
\]

if and only if `T(y) <= C I` for almost every `y`.

**Proof.**  The reverse implication is immediate.  For the forward one, fix
`v` and a compactly supported smooth `eta`, and set

\[
 f_k(y)=k^{-1}\eta(y)\sin(k\langle v,y\rangle).
\]

The leading part of the gradient is
`eta(y) cos(k<v,y>)v`; the remaining part is `O(k^{-1})`.  Letting `k` tend to
infinity in (3.1), by truncation and the Riemann--Lebesgue lemma,

\[
 \int\eta^2 v^T T v\,\rho\leq C|v|^2\int\eta^2\rho.
\]

Lebesgue differentiation, followed by a countable dense set of `v`, proves
the claim.  `square`

Thus the two-step argument

\[
 \operatorname {Var}f
 \leq \mathbb E\langle\tau\nabla f,\nabla f\rangle
 \leq C\mathbb E|\nabla f|^2                         \tag{3.2}
\]

is exactly an `L^infinity` bound on the moment-map Hessian.  It is much
stronger than KLS and is false for measures whose Poincare constants are
uniformly bounded.

By contrast, the assertion

\[
 \operatorname {Var}_\mu f\leq C\mathbb E_\mu|\nabla f|^2 \tag{3.3}
\]

is the Poincare form of KLS itself.  Calling (3.3) a removal lemma without
giving additional content is circular; calling (3.2) a removal lemma is a
strictly stronger and false assertion.

## 4. Spectral truncation and the failure of form-level tail absorption

For `L>0`, use spectral calculus pointwise to write

\[
 \tau=\tau\wedge LI+(\tau-LI)_+.
\]

Then (2.7) gives the exact truncation

\[
 \operatorname {Var}f
 \leq L\mathbb E|\nabla f|^2
 +\mathbb E\langle(\tau-LI)_+\nabla f,\nabla f\rangle. \tag{4.1}
\]

A tempting absorption hypothesis is

\[
 \mathbb E\langle(\tau-LI)_+\nabla f,\nabla f\rangle
 \leq\theta\operatorname {Var}f+B\mathbb E|\nabla f|^2,
 \qquad \theta<1.                                     \tag{4.2}
\]

If true, it would give

\[
 C_P(\mu)\leq {L+B\over1-\theta}.                     \tag{4.3}
\]

For all smooth tests, however, (4.2) again forces
`(tau-LI)_+ <= B I` almost everywhere.  Indeed, apply the oscillatory tests in
Lemma 3.1; their variance is `O(k^{-2})` while their gradient energy has a
nonzero limit.  Thus no distributional tail estimate for `||tau||`, no matter
how strong, can prove (4.2) for arbitrary gradients when `tau` is unbounded.
Gradients can be localized on the bad set and oscillated there.

This leaves two honest possibilities.

* Avoid the weighted upper bound on the bad region and prove a Hardy or
  conditional-variance estimate directly, as happens for a product
  exponential.
* Prove (4.2) only for a low-frequency spectral subspace of the ordinary
  generator.  The oscillatory counterexample then no longer applies.

For the second possibility, let `f` be a normalized first-eigenfunction
approximant, so `Var(f)=1` and `E|nabla f|^2=lambda`.  A universal estimate

\[
 \mathbb E\langle(\tau-LI)_+\nabla f,\nabla f\rangle
 \leq\theta+B\lambda,
 \qquad\theta<1,                                      \tag{4.4}
\]

combined with (4.1) gives `lambda >= (1-theta)/(L+B)`.  This is a precise
noncircular closing lemma, but proving it uniformly would solve KLS.  Neither
thin shell nor slicing presently controls the location and direction of the
gradient of such an eigenfunction relative to the high eigenspaces of `tau`.

## 5. Mandatory model audit

| Isotropic target | Canonical moment kernel | Operator behavior | Test outcome |
|---|---|---|---|
| Gaussian | `I` | identically one | sharp closure |
| cube | `diag((3-X_i^2)/2)` | at most `3/2` | direct closure |
| shifted exponential product | `diag(Y_i)` | mean top eigenvalue `H_d` | `max_i Y_i` loses `log d` |
| regular simplex | `(d+2)(diag(p)-pp^T)` on the tangent space | mean top eigenvalue `Theta(log d)` | `max_i p_i`, rescaled to be 1-Lipschitz, loses `log d` |
| radial exponential | `a(r)P_rad+b(r)P_tan` | radial eigenvalue unbounded | localized radial oscillations rule out form comparison |

### 5.1 Gaussian: exact equality

For the standard Gaussian, take

\[
 \phi(x)={|x|^2\over2}+{d\over2}\log(2\pi).
\]

Then `nabla phi=x`, `tau=I`, and (2.7) is the sharp Gaussian Poincare
inequality.  All trace, determinant, and truncation statements are equalities
in their best possible form.

### 5.2 Isotropic cube: a bounded moment Hessian

First work in one dimension with the uniform probability on `(-a,a)`.  Set

\[
 \phi_a(t)=2\log\cosh(at/2)+\log(4/a).                 \tag{5.1}
\]

Then

\[
 e^{-\phi_a(t)}={a\over4}\operatorname {sech}^2(at/2),
 \quad
 \phi_a'(t)=a\tanh(at/2),
 \quad
 \phi_a''(t)={a^2\over2}\operatorname {sech}^2(at/2). \tag{5.2}
\]

Changing variables `y=phi_a'(t)` shows directly that its density is
`1/(2a)` on `(-a,a)`, and

\[
 \tau(y)={a^2-y^2\over2}.                              \tag{5.3}
\]

For the isotropic cube take `a=sqrt(3)` and tensorize:

\[
 \tau(x)=\operatorname {diag}\left({3-x_1^2\over2},\ldots,
                                    {3-x_d^2\over2}\right),
 \qquad 0\preceq\tau\preceq {3\over2}I.              \tag{5.4}
\]

Thus the moment-map argument itself gives `C_P <= 3/2` for the cube.  The
target density jumps at the boundary, but (5.1)--(5.3) prove the identities
without differentiating that density.  Equivalently, one may approximate the
indicator potential of the cube by smooth convex potentials, isotropize, and
pass to the limit in (2.7).  Any general proof using a smooth target must keep
its estimates uniform in precisely this approximation.

For reference,

\[
 \mathbb E\log\tau_1=\log(3/2)+\int_0^1\log(1-u^2)du
 =\log 6-2.                                            \tag{5.5}
\]

So even this favorable bounded model has an exponentially small typical
determinant.  Determinant size is not the right one-sided datum for removing
the weight.

### 5.3 Product shifted exponentials: the definitive scalar-summary no-go

Let `Y_i` be independent `Exp(1)` variables and put `X_i=Y_i-1`.  Then `X` is
centered and isotropic.  In one dimension the exact moment potential is

\[
 \phi(t)=e^t-t,
 \quad e^{-\phi(t)}=e^{-e^t+t},
 \quad \phi'(t)=e^t-1,
 \quad \phi''(t)=e^t.                                 \tag{5.6}
\]

The substitution `y=e^t` proves normalization and the exponential
pushforward.  Tensorization gives the canonical moment kernel

\[
 \boxed{\ \tau(X)=\operatorname {diag}(Y_1,\ldots,Y_d).\ } \tag{5.7}
\]

It has exceptionally strong scalar control:

\[
\begin{aligned}
 &\mathbb E\tau=I,\\
 &\operatorname {Tr}\tau=\sum_iY_i\sim\operatorname {Gamma}(d,1),
   \qquad \operatorname {Var}(\operatorname {Tr}\tau)=d,\\
 &\mathbb E\log\det\tau=d\,\mathbb E\log Y_1=-\gamma d,
   \qquad \operatorname {Var}(\log\det\tau)={\pi^2\over6}d. \tag{5.8}
\end{aligned}
\]

Nevertheless,

\[
 \|\tau\|_{op}=M_d:=\max_iY_i,
 \quad \mathbb E M_d=H_d,
 \quad \operatorname {Var}M_d=H_d^{(2)}\leq{\pi^2\over6}. \tag{5.9}
\]

The equalities follow from the independent exponential spacings of the order
statistics.  In particular `E||tau||_op ~ log d`, while
`lambda_min(tau)=min_iY_i` has mean `1/d`.  Both the arithmetic data in the
trace and the geometric data in the determinant coexist with severe spectral
anisotropy.

This target already satisfies the strongest relevant target-side scalar
inputs.  Its density has essential supremum one in isotropic position, so its
isotropic constant is one.  Also

\[
 \mathbb E(|X|^2-d)^2
 =d\operatorname {Var}((Y_1-1)^2)=8d,                 \tag{5.10}
\]

because `E(Y-1)^4=9`; hence

\[
 \mathbb E(|X|-\sqrt d)^2\leq8.                       \tag{5.11}
\]

Thus slicing and dimension-free thin shell cannot, by themselves, turn
(5.8) into an operator estimate.

There is also an exact gradient-correlation test.  Let

\[
 f(X)=\max_iY_i.
\]

Away from ties, `|nabla f|=1` and its gradient is the coordinate vector of the
largest entry.  Hence

\[
 \mathbb E|\nabla f|^2=1,
 \quad \operatorname {Var}f=H_d^{(2)}\leq{\pi^2\over6},
 \quad \mathbb E\langle\tau\nabla f,\nabla f\rangle=H_d. \tag{5.12}
\]

The weighted Brascamp--Lieb bound loses `log d` on a natural 1-Lipschitz
function even though the actual variance is bounded.  For fixed `L`, its
tail term in (4.1) is

\[
 \mathbb E(M_d-L)_+,
\]

which is of order `log d-L` until `L` reaches order `log d`.  Choosing
`L ~ log d` merely moves the logarithm into the bounded part of (4.1).

For a purely local no-go, fix a nonzero `eta in C_c^infty((0,1))` and set

\[
 f_{R,k}(X)=k^{-1}\eta(Y_1-R)\sin(kY_1).               \tag{5.13}
\]

As `k` tends to infinity, its weighted-to-unweighted energy ratio is at least
`R+o(1)`, while its variance divided by its unweighted energy tends to zero.
This explicitly disproves both a global form comparison and form-level tail
absorption.

The failure is in the certificate, not in KLS for this model.  The
one-dimensional exponential has Cheeger constant one: if `p(y)=e^{-y}` and
`F(y)=1-e^{-y}`, then

\[
 {p(y)\over\min(F(y),1-F(y))}\geq1.
\]

The one-dimensional Cheeger inequality gives `C_P <= 4`, and tensorization
gives

\[
 C_P\left(\bigotimes_{i=1}^d(\operatorname {Exp}(1)-1)\right)\leq4. \tag{5.14}
\]

This is the prototype of valid Hardy compensation: condition on all other
coordinates, control the one-dimensional variance directly, and then
tensorize.  It never attempts to upper-bound the form in (5.7).

### 5.4 Isotropic simplex: the same logarithm without product target geometry

Let `m=d+1`, let

\[
 H=\{z\in\mathbb R^m:\sum_i z_i=0\},
 \qquad u=m^{-1}(1,\ldots,1),
\]

and let `p=(p_i)` be uniform on the probability simplex, equivalently
`p ~ Dirichlet(1,...,1)`.  The isotropic regular simplex in `H` is

\[
 X=\sqrt{m(m+1)}(p-u).                                 \tag{5.15}
\]

Put `c=sqrt((m+1)/m)` and, up to its normalizing additive constant, define

\[
 \phi(z)=m\log\left(\sum_{i=1}^m e^{c z_i}\right),
 \qquad z\in H.                                       \tag{5.16}
\]

Writing `p_i=e^{cz_i}/sum_j e^{cz_j}`, direct differentiation gives

\[
 \nabla_H\phi=\sqrt{m(m+1)}(p-u),
 \qquad
 \nabla_H^2\phi=(m+1)(\operatorname {diag}p-pp^T)|_H. \tag{5.17}
\]

The softmax Jacobian shows that `e^{-phi}dz` pushes to constant Lebesgue
density on the simplex.  Therefore

\[
 \boxed{\ \tau(X)=(m+1)(\operatorname {diag}p-pp^T)|_H.\ } \tag{5.18}
\]

Since `E pp^T` is scalar on `H`, one checks directly that `E tau=I_H`.
Moreover,

\[
 \operatorname {Tr}_H\tau=(m+1)\left(1-\sum_i p_i^2\right),
 \qquad \mathbb E\operatorname {Tr}_H\tau=m-1=d.     \tag{5.19}
\]

The matrix-tree identity gives

\[
 \det_H(\operatorname {diag}p-pp^T)=m\prod_{i=1}^m p_i,
\]

and therefore

\[
 \mathbb E\log\det_H\tau
 =(m-1)\log(m+1)+\log m-mH_{m-1}
 =-\gamma d+O(1).                                     \tag{5.20}
\]

Let `M=max_i p_i`.  Since

\[
 \lambda_{max}(\operatorname {diag}p-pp^T)\leq M
\]

and, for an index attaining `M`, the unit vector proportional to `e_i-u`
has Rayleigh quotient

\[
 {m\over m-1}M(1-M),                                  \tag{5.21}
\]

the expected top eigenvalue in (5.18) is `Theta(log m)`.  Indeed
`E M=H_m/m`, and the contribution of `M>1/2` is exponentially small.

The gradient correlation can again be made exact.  Define

\[
 F(X)=m\sqrt{{m+1\over m-1}}\max_i p_i.               \tag{5.22}
\]

Then `|nabla_H F|=1` away from ties, while

\[
 \mathbb E\langle\tau\nabla F,\nabla F\rangle
 ={m(m+1)\over m-1}\mathbb E[M(1-M)]
 =H_m+O\left({H_m^2+1\over m}\right).                 \tag{5.23}
\]

Its variance is universally bounded.  To see this without an asymptotic
theorem, write `p_i=E_i/S`, where the `E_i` are independent exponentials and
`S=sum E_i`.  The proportions `p` are independent of `S`.  If
`A=max_i E_i=SM`, then

\[
 \operatorname {Var}A=H_m^{(2)}\leq{\pi^2\over6}
 =m(m+1)\operatorname {Var}M+m(\mathbb EM)^2.
\]

It follows that

\[
 \operatorname {Var}F
 ={m+1\over m-1}\operatorname {Var}(mM)
 \leq {m+1\over m-1}{\pi^2\over6}.                   \tag{5.24}
\]

Thus the simplex reproduces the exact pattern of (5.12): unit unweighted
energy, bounded variance, and logarithmic moment-map energy.  The obstruction
is not merely tensor-product geometry.

### 5.5 Isotropic radial exponential: an unbounded canonical radial eigenvalue

Let

\[
 d\mu_n(y)=Z_n^{-1}e^{-\alpha|y|}dy,
 \qquad \alpha=\sqrt{n+1}.                            \tag{5.25}
\]

Its radius `S=|Y|` is `Gamma(n,alpha)`, so

\[
 \mathbb ES^2={n(n+1)\over\alpha^2}=n,
 \qquad \operatorname {Var}S={n\over n+1}<1.          \tag{5.26}
\]

Hence the measure is isotropic and has an exact dimension-free thin shell.
Stirling's formula also gives a universal isotropic constant.

Rotational invariance and uniqueness of the moment potential allow it to be
chosen radial, `phi(x)=h(r)`, `r=|x|`.  Put

\[
 s=h'(r),\qquad a(s)=h''(r),\qquad b(s)={h'(r)\over r}={s\over r}. \tag{5.27}
\]

At `y=s theta`, the canonical kernel is exactly

\[
 \tau(y)=a(s)\,\theta\theta^T
          +b(s)(I-\theta\theta^T).                    \tag{5.28}
\]

The radial Monge--Ampere equation is

\[
 a(s)b(s)^{n-1}=C_n\exp[-h(r)+\alpha s].              \tag{5.29}
\]

This implicit formula is enough for a useful no-go.

**Lemma 5.1.**  The radial eigenvalue `a(s)` of the canonical moment kernel is
unbounded.

**Proof.**  The target radius is unbounded, so `s=h'(r)` is unbounded.  Suppose
instead that `a<=A` eventually.  If the radial domain of `h` has a finite
endpoint, boundedness of `a=ds/dr` already makes `s` bounded, a contradiction.
The radial domain must therefore be unbounded.  Once `s>=2 alpha A`,
differentiating the logarithm of (5.29) with respect to `r` gives

\[
 {d\over dr}\log(ab^{n-1})=-s+\alpha a\leq-s/2.
\]

Thus `ab^{n-1} <= C e^{-cr}` from some point onward.  Since `s` is then
bounded below by a positive constant,

\[
 a=(ab^{n-1})(r/s)^{n-1}\leq C r^{n-1}e^{-cr}.
\]

The right side is integrable, which makes
`s(r)=s(R)+int_R^r a(u)du` bounded, a contradiction.  `square`

Consequently direct form comparison fails in this rotationally invariant
model too.  If `a(s_0)` is large, continuity gives a radial interval on which
it remains large.  The explicit radial oscillatory test

\[
 f_k(y)=k^{-1}\eta(|y|)\sin(k|y|),                    \tag{5.30}
\]

with `eta` supported in that interval makes the weighted-to-unweighted energy
ratio arbitrarily large while its variance-to-energy ratio tends to zero.

For comparison, this measure has another particularly simple Stein kernel,
not asserted to carry the moment-map Brascamp--Lieb inequality:

\[
 \widehat\tau(y)=\left({|y|\over\alpha}+{1\over\alpha^2}\right)I. \tag{5.31}
\]

Indeed, if `rho(y)=Z_n^{-1}e^{-alpha|y|}` and
`w(r)=r/alpha+alpha^{-2}`, then

\[
 \operatorname {div}(\rho wI)
 =\rho(w'-\alpha w){y\over|y|}=-y\rho.
\]

Also `E w=1`.  Thus rotational scalarization can simplify a Stein identity,
but a Stein identity alone is not the weighted Poincare certificate (2.7),
and even this scalar kernel is unbounded.

## 6. What determinant and trace information cannot do

The product-exponential matrix in (5.7) is an exact realized matrix model with

\[
 \mathbb E\tau=I,
 \quad \operatorname {Tr}\tau=d+O_{L^2}(\sqrt d),
 \quad {1\over d}\log\det\tau=-\gamma+O_{L^2}(d^{-1/2}),
 \quad \|\tau\|_{op}\asymp\log d.                    \tag{6.1}
\]

Therefore none of the following implications is valid:

* concentrated trace plus `E tau=I` implies bounded top eigenvalue;
* a dimension-free lower bound on the geometric mean of the eigenvalues
  implies bounded top eigenvalue;
* thin shell of the target controls the top eigenvalue of its moment kernel;
* the slicing bound on the target density controls that eigenvalue.

For fixed `L`, the high-rank statistic is also explicit:

\[
 \mathbb E\operatorname {rank}1_{\{\tau>L\}}=de^{-L},
 \qquad
 \mathbb E\operatorname {Tr}(\tau-LI)_+=de^{-L}.      \tag{6.2}
\]

At each point only a small fraction of directions may be bad, but a nonlinear
gradient can select a bad direction adaptively.  The max test (5.12) does
exactly that.  Averaging projection estimates or using a low expected bad
rank cannot ignore this selection.

The simplex shows the same phenomenon with a non-diagonal kernel.  There the
high direction is the contrast between the largest barycentric coordinate and
the remaining coordinates, and the gradient of the largest-coordinate
function selects precisely that contrast.

## 7. Tail/Hardy compensation: what works and what cannot work

The shifted exponential separates two notions which are easy to conflate.

1. **False form compensation.**  Bound the large part of
   `E<tau nabla f,nabla f>` by ordinary energy and a fraction of the variance.
   Lemma 3.1 and (5.13) rule this out.
2. **Valid variance compensation.**  Return to a conditional variance before
   applying weighted Brascamp--Lieb.  On each one-dimensional exponential
   fiber, use the Cheeger/Hardy inequality directly; then tensorize.  This is
   (5.14).

For a general moment kernel, the second program would require a geometric
decomposition into fibers on which large eigenvalues have a one-dimensional
Hardy interpretation, together with a dimension-free rule for recombining
the conditional variances.  A spectral decomposition of `tau(y)` is not such
a foliation: its eigenspaces vary with `y`, can exchange order, and need not
integrate to global fibers.  Although `tau^{-1}=nabla^2 phi^*` in the smooth
invertible case supplies Hessian structure, Monge--Ampere controls the product
of eigenvalues, not the global geometry of these prospective fibers.

The max tests show what a successful compensation must recognize.  In the
product case, the identity of the largest coordinate changes across the
codimension-one tie sets, and conditional one-dimensional inequalities exploit
all coordinates instead of charging the selected maximum weight.  In the
simplex the analogous chambers meet under the barycentric constraint.  Trace,
determinant, and tail probabilities contain none of this chamber-gluing data.

## 8. Randomized transport or kernel ensembles

There are three different averaging operations.

### 8.1 Symmetries of the canonical moment map give no new kernels

The centered moment potential is unique up to source translation.  If `Q` is
an orthogonal symmetry of `mu`, equivariance gives

\[
 \tau_Q(y)=Q\tau(Q^Ty)Q^T=\tau(y).                    \tag{8.1}
\]

Thus random source translations, target rotations followed by rotation back,
or actual symmetry transformations reproduce the same canonical kernel.  For
the product exponential, coordinate permutations in (8.1) leave
`diag(Y_i)` exactly unchanged as a field.

### 8.2 Conjugating the matrix while freezing the point is invalid

If one independently permutes the diagonal entries of (5.7) while leaving
`X` fixed, the average is

\[
 \bar\tau(X)={\sum_iY_i\over d}I.                     \tag{8.2}
\]

This looks favorable but is neither a Stein kernel nor a weighted Poincare
weight.  For `g(X)=X_1^2/2`,

\[
 \mathbb E[X_1g(X)]={1\over2}\mathbb EX_1^3=1,
 \qquad
 \mathbb E[\bar\tau_{11}\partial_1g]={1\over d}.      \tag{8.3}
\]

And

\[
 \operatorname {Var}(X_1^2/2)=2,
 \qquad
 \mathbb E\langle\bar\tau\nabla g,\nabla g\rangle
 =1+{2\over d}<2\quad(d>2).                           \tag{8.4}
\]

The coupling between the matrix and its spatial argument is indispensable.

### 8.3 A genuine bounded ensemble would be a new KLS certificate

Suppose one constructed weights `tau_omega` such that, for every `omega`,

\[
 \operatorname {Var}_\mu f
 \leq\mathbb E_\mu\langle\tau_\omega\nabla f,\nabla f\rangle,
\]

and their pointwise average obeyed

\[
 \mathbb E_\omega\tau_\omega(y)\preceq CI.            \tag{8.5}
\]

Averaging would prove `C_P(mu)<=C`.  This is a legitimate possible route, not
a logical contradiction.  But positivity and the Stein identity of an
averaged matrix are not enough: each member must carry a variance inequality,
and (8.3)--(8.4) show that naive rotational decorrelation destroys both
properties.  No such ensemble follows from moment-map uniqueness, thin shell,
or slicing.  Establishing (8.5) uniformly would itself be the new central
ingredient.

## 9. Circularity and implication ledger

The precise logical status of the common candidate assumptions is as follows.

* `tau(y) <= C I` implies KLS through (2.7), but it is strictly stronger and
  false for the shifted exponential and radial exponential moment maps.
* `E<tau nabla f,nabla f> <= C E|nabla f|^2` for every `f` is equivalent to
  the same pointwise bound by Lemma 3.1, so it is also false.
* `E tau=I` is exact and tests only constant gradients.  It gives no control
  of gradients selected from the location of the large eigenvalues.
* `E||tau||_op <= C` would be useful only with an additional decorrelation
  argument, and in any event is false for both the product exponential and
  simplex, where it is `Theta(log d)`.
* Bounds on `E log det tau`, `Tr tau`, their concentration, or the expected
  rank of a spectral tail do not close the gradient correlation.  Equations
  (5.8), (5.12), and (6.2) are an explicit simultaneous countermodel.
* A direct universal inequality `Var f <= C E|nabla f|^2` is exactly the
  Poincare/KLS target.  Deriving it by invoking an ordinary semigroup spectral
  gap, a resolvent norm, or a Hardy inequality whose best constant is
  `C_P(mu)` is circular.
* The low-spectrum estimate (4.4) is not formally circular.  It is a sharply
  stated sufficient lemma which avoids all high-frequency no-go tests.
  However, a uniform proof of it is conjecture-strength and requires input
  beyond the currently used scalar consequences of thin shell and slicing.
* A valid randomized family satisfying (8.5) is likewise a genuine sufficient
  certificate.  Canonical moment-map symmetries do not produce such a family,
  and independent matrix randomization is invalid.

## 10. Remaining viable formulation

The moment-map route should not pursue a pointwise matrix bound, an averaged
operator norm, or a form bound for arbitrary gradients.  All three are
decisively excluded by explicit log-concave models.

The narrow viable question is instead the following low-mode correlation
statement.

> Let `mu` be isotropic and log-concave, let `tau` be its moment-map kernel,
> and let `f` lie in a sufficiently low spectral subspace of the ordinary
> Dirichlet form.  Can dimension-free thin shell and slicing, together with
> the Hessian/Monge--Ampere structure of `tau`, prevent `nabla f` from carrying
> most of its weighted energy in the high eigenspaces of `tau`?

Quantitatively, (4.4) is one adequate answer.  The cube and Gaussian pass it
trivially.  The exponential, simplex, and radial high-frequency tests do not
refute the spectral restriction, while they show why that restriction is
essential.  A proof would need new information coupling an ordinary low mode
to the spatial geometry of the moment-map Hessian; trace, determinant, target
thin shell, target density height, and symmetry averaging do not provide that
coupling.
