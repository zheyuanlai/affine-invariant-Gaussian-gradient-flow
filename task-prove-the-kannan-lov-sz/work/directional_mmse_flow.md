# Directional MMSE at one Gaussian unit of noise

## 0. Verdict

Let \(X\sim\mu\) be isotropic and log-concave in \(\mathbb R^n\), let
\(G\sim N(0,I_n)\) be independent, and put

\[
 Y=X+G,\qquad
 A(Y)=\operatorname {Cov}(X\mid Y).
\]

The proposed directional estimate is

\[
 \boxed{
 \inf_{|b|=1}\mathbb E\,b^TA(Y)b\ge \kappa
 }
 \tag{DMMSE}
\]

for a numerical \(\kappa>0\), uniformly in dimension and in \(\mu\).

I found neither an unconditional proof nor a log-concave counterexample.
The exact audit is as follows.

1. KLS implies (DMMSE) immediately.  In fact

   \[
   \mathbb E\,b^TA(Y)b
   \ge {1\over C_P(\mu)+2}.
   \tag{0.1}
   \]

   Thus a sequence violating (DMMSE) would also be a sequence with
   \(C_P(\mu)\to\infty\), hence a counterexample to KLS.  The converse
   implication, from (DMMSE) to the full KLS conjecture, is not established;
   the statement is a directional Gaussian-channel consequence of KLS, not
   presently a proved reformulation of all of KLS.

2. Entropy gives the exact **trace** estimate

   \[
   \mathbb E\operatorname {Tr}A(Y)
   \ge {n\over2}N(X),
   \qquad
   N(X):={1\over2\pi e}e^{2h(X)/n}.
   \tag{0.2}
   \]

   For an isotropic log-concave density with isotropic constant
   \(L_\infty=\|f\|_\infty^{1/n}\), this yields

   \[
   \mathbb E\operatorname {Tr}A(Y)
   \ge {n\over4\pi eL_\infty^2}.
   \]

   A universal linear trace floor from this route requires a universal
   bound on \(L_\infty\), namely the slicing/hyperplane conjecture.  The
   classical unconditional estimate \(L_\infty\lesssim n^{1/4}\) gives only
   a \(c\sqrt n\) trace floor.  Even a conjectural \(cn\) trace floor would
   not give (DMMSE): it controls the average posterior error over an
   orthonormal basis, whereas (DMMSE) excludes one exceptional recoverable
   direction.  Fisher information expresses this gap exactly:

   \[
   \mathbb E A(Y)=I-J(\mu*\gamma).
   \tag{0.3}
   \]

   Trace entropy controls \(\operatorname {Tr}J\), while (DMMSE) asks for the
   operator inequality \(J(\mu*\gamma)\preceq(1-\kappa)I\).

3. Products and Gaussians satisfy (DMMSE).  Irreducible orthogonal symmetry
   reduces the directional question to the trace question; when combined
   with a universal entropy-power bound, this covers regular simplices,
   cubes, cross-polytopes, and Euclidean balls.  Irreducibility alone should
   not be mistaken for a proof of slicing for every invariant density.
   Explicit skew-cone, half-ball, variable-width wedge, and thin-tube stress
   tests do not produce a vanishing direction.  In each many-coordinate
   encoding attempt, log-concavity concentrates the longitudinal variable
   at exactly the scale that keeps the aggregate transverse information
   bounded.

4. If (DMMSE) were available, it would close the near-linear eigenfunction
   branch **directly**, without a final-time alignment theorem.  If

   \[
   f=a\cdot x+r,\qquad \|r\|_2=\delta,\qquad
   |a|=\sqrt{1-\delta^2},\qquad b={a\over|a|},
   \]

   and \(c_1=\operatorname {Cov}_{\mu_1}(X,f)\), then

   \[
   c_1=A_1a+e_1,qquad
   e_1=\operatorname {Cov}_{\mu_1}(X,r),
   \]

   and the posterior Poincare inequality gives

   \[
   \boxed{
   \sqrt\lambda
   \ge \|c_1\|_{L^2}
   \ge \kappa\sqrt{1-\delta^2}-\delta.}
   \tag{0.4}
   \]

   Thus, for example, \(\delta\le\kappa/4\) implies
   \(\lambda\ge\kappa^2/4\).

5. The floor does **not** imply the proposed final alignment estimate from
   the terminal identities involving

   \[
   J_t=A_tM_t,qquad c_1=A_1a+e_1.
   \]

   Even when \(e_1=0\), that alignment quotient contains a negative weight
   \((b^TA_1^2b)^{-1}\).  An average lower bound on \(b^TA_1b\) gives no
   control of rare events on which this denominator is small and \(M_1b\)
   is large.  Section 7 gives a terminal matrix model satisfying the MMSE
   floor in every direction, \(\mathbb EJ_1=I\), \(\mathbb EM_1=2I\), and
   the trace scale \(\mathbb E\operatorname {Tr}(M_1^TM_1)=O(n)\), while
   the selected alignment is of order \(n\).  The model is an algebraic
   no-go for a conversion based only on these endpoint facts; it is not
   asserted to be an actual stochastic-localization path.

The most useful conclusion is therefore conditional but sharp:

> A universal directional MMSE floor would solve the near-linear branch by
> the elementary endpoint estimate (0.4).  Entropy supplies only a trace
> estimate at the entropy-power/isotropic-constant scale; making even that
> trace estimate linear invokes slicing, and upgrading a linear trace
> statement to a direction chosen by the eigenfunction is a further
> unresolved operator step.

---

## 1. Exact Gaussian-channel identities

Let \(q\) be the density of \(Y\), and define

\[
 m(y)=\mathbb E[X\mid Y=y],\qquad
 A(y)=\operatorname {Cov}(X\mid Y=y).
 \tag{1.1}
\]

The posterior is

\[
 d\mu_y(x)
 ={\exp(-|x-y|^2/2)\over
   \int\exp(-|z-y|^2/2)d\mu(z)}d\mu(x).
 \tag{1.2}
\]

It is \(1\)-strongly log-concave on its convex support, so Brascamp--Lieb
gives

\[
 0\preceq A(y)\preceq I.                         \tag{1.3}
\]

Differentiation of the exponential family gives the matrix Tweedie
identities

\[
 Dm(y)=A(y),\qquad
 \nabla\log q(y)=m(y)-y,\qquad
 D^2\log q(y)=A(y)-I.                              \tag{1.4}
\]

Let \(J(q)\) denote the Fisher-information matrix

\[
 J(q)=\mathbb E[\nabla\log q(Y)\nabla\log q(Y)^T].
 \tag{1.5}
\]

Since \(m(Y)-Y=-\mathbb E[G\mid Y]\) and
\(\operatorname {Cov}(G\mid Y)=A(Y)\), conditional variance gives

\[
 \boxed{\mathbb E A(Y)=I-J(q).}                    \tag{1.6}
\]

For a unit vector \(b\), set

\[
 \varepsilon_b
 :=\mathbb E\operatorname {Var}(b\cdot X\mid Y)
 =\mathbb E\,b^TA(Y)b.                             \tag{1.7}
\]

Then (1.6) says

\[
 \boxed{
 \varepsilon_b=1-\mathbb E(\partial_b\log q(Y))^2.}
 \tag{1.8}
\]

Thus (DMMSE) is precisely a uniform spectral gap below the Gaussian-noise
ceiling for the directional Fisher information of \(q=\mu*\gamma\).

There is also a useful posterior-mean form.  Put

\[
 g_b(y)=b\cdot m(y)=\mathbb E[b\cdot X\mid Y=y].
 \tag{1.9}
\]

Total variance and (1.4) give

\[
 \operatorname {Var}_q(g_b)=1-\varepsilon_b,
 \qquad
 \nabla g_b(y)=A(y)b.                              \tag{1.10}
\]

Because \(0\preceq A\preceq I\),

\[
 \mathbb E|\nabla g_b|^2
 =\mathbb E\,b^TA^2b
 \le\mathbb E\,b^TAb
 =\varepsilon_b.                                  \tag{1.11}
\]

These identities are valid by smooth approximation for arbitrary
full-dimensional log-concave measures, including uniform measures on
convex bodies.

### 1.1 A pointwise lower bound is false

Only the averaged statement can be plausible.  If \(X\) is uniform on a
bounded interval, then for observations \(y\) far beyond an endpoint the
posterior concentrates in a boundary layer whose variance tends to zero.
Thus

\[
 A(y)\succeq \kappa I\quad\hbox{for every }y
 \tag{1.12}
\]

is already false in dimension one.  Any proof must use the output average
in (1.7).

---

## 2. The exact Poincare/KLS implication

Apply the Poincare inequality of \(q\) to \(g_b\).  Equations
(1.10)--(1.11) give

\[
 1-\varepsilon_b
 \le C_P(q)\mathbb E|\nabla g_b|^2
 \le C_P(q)\varepsilon_b.
 \tag{2.1}
\]

Therefore

\[
 \boxed{
 \varepsilon_b\ge {1\over 1+C_P(q)}.}
 \tag{2.2}
\]

Poincare constants subadd under convolution:

\[
 C_P(\mu*\gamma)\le C_P(\mu)+C_P(\gamma)
 =C_P(\mu)+1.                                      \tag{2.3}
\]

For completeness, condition a test function \(h(X+G)\) first on \(G\),
apply Poincare for \(\mu\), then condition on \(X\), apply Gaussian
Poincare, and use Jensen on the two conditional gradients.  Combining
(2.2)--(2.3) gives

\[
 \boxed{
 \varepsilon_b\ge {1\over C_P(\mu)+2}.}
 \tag{2.4}
\]

This has two immediate logical consequences.

* **[KLS INPUT]** Substituting a universal bound
  \(C_P(\mu)\le C\) in (2.4) is exactly the Poincare formulation of the KLS
  conjecture.  This proves (DMMSE) conditionally on KLS.

* If \(\varepsilon_{b_n}\to0\) for a sequence \((\mu_n,b_n)\), then

  \[
  C_P(\mu_n*\gamma)\ge {1-\varepsilon_{b_n}\over
                              \varepsilon_{b_n}}
  \quad\hbox{and hence}\quad
  C_P(\mu_n)\to\infty.                             \tag{2.5}
  \]

  Consequently, an explicit counterexample to (DMMSE) would refute KLS.

The reverse implication is not contained in these formulas.  A lower
bound for the posterior error of every **linear** statistic does not
immediately give Poincare or Cheeger control for arbitrary functions or
sets.  Accordingly, (DMMSE) should be labelled KLS-relevant and implied by
KLS, but not claimed to be KLS-equivalent without a new reduction.

It is equally circular to assert directly that

\[
 J(\mu*\gamma)\preceq(1-\kappa)I:                 \tag{2.6}
\]

by (1.6), this is just (DMMSE) in Fisher notation.

---

## 3. What entropy proves, and where it stops

Define the entropy power per coordinate by

\[
 N(X):={1\over2\pi e}\exp\left({2h(X)\over n}\right).
 \tag{3.1}
\]

It is important not to insert a universal lower bound for \(N(X)\) here.
For general isotropic log-concave vectors, a dimension-free lower bound on
\(N(X)\) is equivalent up to numerical constants to the
slicing/hyperplane conjecture.  The entropy calculation below is exact,
but it retains \(N(X)\).

Let

\[
 D:=\mathbb E\operatorname {Tr}A(Y)
 =\mathbb E|X-m(Y)|^2.                              \tag{3.2}
\]

Since \(\operatorname {Cov}(Y)=2I\), Gaussian maximal entropy gives

\[
 h(Y)\le {n\over2}\log(4\pi e).                   \tag{3.3}
\]

The additive channel identity and (3.1)--(3.3) yield

\[
\begin{aligned}
 h(X\mid Y)
 &=h(X)+h(G)-h(Y)\\
 &\ge {n\over2}\log(\pi e N(X)).                 \tag{3.4}
\end{aligned}
\]

On the other hand, with \(E=X-m(Y)\), conditional translation followed by
Gaussian maximal entropy and determinant--trace AM--GM gives

\[
\begin{aligned}
 h(X\mid Y)
 &=h(E\mid Y)\\
 &\le h(E)\\
 &\le {n\over2}\log\left(2\pi e{D\over n}\right).
                                                               \tag{3.5}
\end{aligned}
\]

Comparison gives the exact entropy trace estimate

\[
 \boxed{
 \mathbb E\operatorname {Tr}A(Y)
 \ge {N(X)\over2}n.}
 \tag{3.6}
\]

This rate--distortion argument is unconditional, but its numerical scale
is the input entropy power.

If \(f\) is the density of an isotropic log-concave \(X\), define

\[
 L_\infty(X):=\|f\|_\infty^{1/n}.                  \tag{3.7}
\]

Since \(h(X)\ge-\log\|f\|_\infty\),

\[
 N(X)\ge {1\over2\pi eL_\infty(X)^2}.
 \tag{3.8}
\]

Consequently,

\[
 \boxed{
 \mathbb E\operatorname {Tr}A(Y)
 \ge {n\over4\pi eL_\infty(X)^2}.}
 \tag{3.9}
\]

A universal bound \(L_\infty(X)\le C\) is the slicing/hyperplane
conjecture in functional form.  Conversely, for log-concave densities the
standard comparison

\[
 -\log\|f\|_\infty
 \le h(X)
 \le-\log\|f\|_\infty+n                           \tag{3.10}
\]

shows that a dimension-free entropy-power lower bound has the same
dimension-free content.  Thus a claim that reverse entropy gives
\(\mathbb E\operatorname {Tr}A\ge cn\) for every isotropic log-concave law
would already be importing slicing.  The classical unconditional bound
\(L_\infty\le Cn^{1/4}\), for example, gives the weaker but valid estimate

\[
 \mathbb E\operatorname {Tr}A(Y)\ge c\sqrt n.      \tag{3.11}
\]

Marsiglietti--Kostina give explicit dimension-dependent vector entropy
bounds and dimension-free constants under additional symmetry such as
unconditionality; their result does not supply a universal entropy-power
lower bound for arbitrary isotropic log-concave vectors.

### 3.1 The operator gap

Equation (3.6) says only

\[
 \operatorname {Tr}\mathbb EA\ge c n.
\]

Even if slicing supplied a linear trace floor, the algebraic matrix

\[
 \overline A=\operatorname {diag}(\epsilon,1,\ldots,1)
 \tag{3.12}
\]

shows why this cannot by itself imply
\(\overline A\succeq c'I\).  Isotropy fixes the **prior** covariance but
does not make the averaged posterior covariance rotationally invariant.
In Fisher form, (3.6) controls \(\operatorname {Tr}J(q)\), whereas
(DMMSE) controls \(\|J(q)\|_{\mathrm{op}}\).

* **[UNCONDITIONAL ENTROPY]** The exact conclusion (3.6), with the factor
  \(N(X)\), is unconditional.

* **[SLICING INPUT]** Replacing \(N(X)\) by a numerical constant, or
  equivalently bounding \(L_\infty(X)\) universally in (3.9), invokes the
  slicing/hyperplane conjecture.  KLS implies this input, but the two
  conjectures should not be identified.

* **[OPEN OPERATOR STEP]** Even after a linear trace bound is available,
  upgrading it to a lower bound in a prescribed direction is a separate
  missing step.  Doing it by inserting a dimension-free Poincare constant
  returns to the KLS input in Section 2.

No determinant estimate helps here without directional control: a single
posterior eigenvalue may be very small while the trace and entropy power
remain of order \(n\).

---

## 4. Classes where the directional floor follows

### 4.1 One dimension and products

In one dimension an isotropic log-concave law has a universal Poincare
constant.  For example, the standard estimate \(C_P\le12\operatorname
{Var}\) and (2.4) give

\[
 \mathbb E\operatorname {Var}(X\mid X+G)\ge {1\over14}.
 \tag{4.1}
\]

If \(\mu=\mu_1\otimes\cdots\otimes\mu_n\) is an isotropic product, the
posterior factorizes coordinatewise.  Hence for \(b=(b_i)\),

\[
 \mathbb E\operatorname {Var}(b\cdot X\mid Y)
 =\sum_i b_i^2\mathbb E\operatorname {Var}(X_i\mid X_i+G_i)
 \ge {1\over14}.                                   \tag{4.2}
\]

Thus high dimension by itself does not create a coding counterexample;
dependence between the selected direction and many transverse coordinates
is essential.

### 4.2 Gaussian input

For \(X\sim N(0,I)\),

\[
 A(y)={1\over2}I,
 \qquad
 \varepsilon_b={1\over2}.                         \tag{4.3}
\]

### 4.3 Irreducible orthogonal symmetry

Suppose an orthogonal group \(\mathcal G\) preserves \(\mu\) and acts
irreducibly on \(\mathbb R^n\).  Channel equivariance implies that
\(\mathbb EA(Y)\) commutes with every element of \(\mathcal G\).  Since it
is self-adjoint, irreducibility forces

\[
 \mathbb EA(Y)=\alpha I.                           \tag{4.4}
\]

The exact trace estimate (3.6) then gives

\[
 \alpha\ge {N(X)\over2}
 \ge {1\over4\pi eL_\infty(X)^2}.                 \tag{4.5}
\]

Thus symmetry removes the **trace-to-direction** obstruction, but it does
not by itself prove a universal entropy-power bound for every invariant
density.  For the standard canonical bodies--Euclidean balls, cubes,
cross-polytopes, and regular simplices in their natural isotropic
positions--the isotropic constants are uniformly bounded, so (4.5) does
give (DMMSE) with a numerical constant.  Cubes and cross-polytopes also
fall under unconditional entropy estimates, and cubes already follow from
the product argument.

This symmetry argument is useful diagnostically.  A counterexample must be
both highly dependent and strongly anisotropic at the level of posterior
recoverability, even though its prior covariance is the identity.

---

## 5. Skew cones and variable-width bodies

The most plausible geometric counterexample would encode one longitudinal
coordinate into the aggregate shape of many transverse coordinates.  The
following calculations show the compensation repeatedly imposed by
log-concavity.

### 5.1 A cone or simplex height

For a cone with homothetic \((n-1)\)-dimensional slices, write the
normalized height from the base as \(U\in[0,1]\).  Uniform volume has

\[
 \rho_U(u)\propto(1-u)^{n-1},
 \qquad U\sim\operatorname {Beta}(1,n).            \tag{5.1}
\]

Thus

\[
 \mathbb EU={1\over n+1},\qquad
 \operatorname {Var}(U)\asymp {1\over n^2}.       \tag{5.2}
\]

The standardized longitudinal variable is therefore of the form
\(T\simeq nU-1\).  Conditional transverse lengths are multiplied by
\(1-U\), so an order-one change in \(T\) changes each transverse variance
by only \(O(n^{-1})\).  An aggregate quadratic statistic of \(n-1\) noisy
transverse observations fluctuates on relative scale \(n^{-1/2}\), which
is larger than this \(n^{-1}\) modulation.  Equivalently, the local Fisher
information from \(n\) coordinates, each with slope \(O(n^{-1})\), is only
\(O(n\cdot n^{-2})=O(n^{-1})\).

For a regular simplex this calculation is exact in barycentric
coordinates.  If \((U_0,\ldots,U_n)\) is uniform on the probability
simplex, then \(U_0\sim\operatorname {Beta}(1,n)\); conditional on \(U_0=u\),
the opposite slice is a copy of the \((n-1)\)-simplex scaled by \(1-u\).
The full symmetry argument in Section 4.3 already gives the rigorous
directional floor.  The height calculation explains why the tempting
``estimate height from transverse radius'' mechanism does not beat the
unit longitudinal Gaussian noise.

The same boundary-layer scaling applies to a highly skew Euclidean cone:
making the cone taller changes the affine normalization of the height, but
after isotropization its typical height fluctuation still corresponds to a
relative transverse-width change of order \(1/n\).

### 5.2 A paired-slope wedge: a quantitative compensation

Consider the convex body

\[
 K_a=\left\{(t,x,y):
 |x_i|\le1+a_it,\quad |y_i|\le1-a_it,
 \quad 1\le i\le m\right\},                        \tag{5.3}
\]

on the interval where all widths are positive.  The longitudinal marginal
is

\[
 \rho(t)\propto\prod_{i=1}^m(1-a_i^2t^2).          \tag{5.4}
\]

Its potential satisfies

\[
\begin{aligned}
 {d^2\over dt^2}\left[-\log\rho(t)\right]
 &=\sum_i{2a_i^2(1+a_i^2t^2)\over(1-a_i^2t^2)^2}\\
 &\ge2\sum_i a_i^2.                                \tag{5.5}
\end{aligned}
\]

If \(v=\operatorname {Var}(t)\), Brascamp--Lieb gives

\[
 v\sum_i a_i^2\le {1\over2}.                      \tag{5.6}
\]

Now standardize \(T=t/\sqrt v\) and each transverse coordinate.  The
linearized conditional-variance slope of the normalized \(x_i\) coordinate
is

\[
 \ell_i={2a_i\sqrt v\over1+a_i^2v},                \tag{5.7}
\]

and the \(y_i\) coordinate has the opposite slope.  Consequently

\[
 \sum_{i=1}^m2\ell_i^2
 \le8v\sum_i a_i^2
 \le4.                                             \tag{5.8}
\]

Thus increasing the number of transverse channels cannot make their total
local variance-modulation signal diverge: the same slice-volume factor
which creates the channels makes the longitudinal marginal more curved.
Adding unit Gaussian noise cannot increase the information available from
the transverse sample.  Estimate (5.8) rigorously establishes the
concentration-versus-slope compensation at the linearized level.  It is not
a complete MMSE lower bound for this body, because a full proof would also
need uniform control of nonlinear likelihood and tail effects, but it
defeats the naive claim that merely adding more width channels makes their
standardized first-order signal diverge.

### 5.3 A curved radial slice

For a half-ball or a body of revolution, a distinguished coordinate \(s\)
has slice volume proportional to a high power of the radius.  In the unit
ball,

\[
 \rho(s)\propto(1-s^2)^{(n-1)/2},
 \qquad s=O(n^{-1/2}).                              \tag{5.9}
\]

After longitudinal standardization \(T\simeq\sqrt n\,s\), the conditional
transverse radius is

\[
 \sqrt{1-s^2}=1-{T^2\over2n}+O(n^{-2}T^4).         \tag{5.10}
\]

Again, its modulation is only \(O(n^{-1})\) across an order-one change of
the standardized height.  The \(n\) transverse coordinates do not supply
diverging information.  A one-sided truncation destroys full irreducible
symmetry but leaves this compensation mechanism intact.

---

## 6. Near-degenerate tubes and nonlinear encodings

### 6.1 Thin support around a curve

Suppose a convex body were contained in an \(\epsilon\)-tube around a
macroscopic graph

\[
 \{(t,\Gamma(t)):t\in I\}.                         \tag{6.1}
\]

Take two points over \(t_0,t_1\).  Convexity puts their entire chord in the
body.  If \(\epsilon\to0\) while \(I\) stays macroscopic, this forces

\[
 \Gamma((1-s)t_0+st_1)
 =(1-s)\Gamma(t_0)+s\Gamma(t_1)                    \tag{6.2}
\]

for every \(s\in[0,1]\).  Hence the limiting graph is affine.  A genuinely
curved macroscopic tube must have thickness comparable to its chordal
deviation; it cannot provide arbitrarily many nearly noiseless nonlinear
measurements of \(t\) while remaining convex.

For a log-concave density the same obstruction appears through convex
sublevel sets of its potential: a very narrow valley over a long interval
cannot track a nonlinear curve without its convex hull filling in the
chords.

### 6.2 Affine tubes

An affine tube can be straightened by a linear map.  Isotropic
normalization then expands its thin directions to unit variance and removes
the linear regression of the transverse coordinates on the longitudinal
one.  Fixed-width affine tubes reduce to product-like tests, for which
Section 4.1 gives a constant floor.  Variable-width affine tubes return to
the slice-volume compensation in Section 5.2.

### 6.3 Perspective models

A smooth analogue makes the same phenomenon transparent.  The function

\[
 (t,z)\longmapsto {|z|^2\over2t},\qquad t>0,        \tag{6.3}
\]

is jointly convex.  Conditional on \(t\), the \(m\) transverse coordinates
behave like Gaussians of variance \(t\), but integration in \(z\) produces
the longitudinal factor \(t^{m/2}\).  Its logarithmic curvature is of
order \(m/t^2\), so the standardized fluctuation of \(t\) is only
\(O(m^{-1/2})\) relative to its mean.  Each transverse variance then has
slope \(O(m^{-1/2})\), and the aggregate information from \(m\) coordinates
stays of constant order rather than diverging.

These tests do not prove the general theorem, but no cone, skew slice, or
near-degenerate curved/tube construction found here evades the same
concentration-versus-encoding balance.

---

## 7. Consequences for the localization flow

At time one in stochastic localization, the random tilt parameter has the
same channel law as \(Y=X+G\), and

\[
 A_1=\operatorname {Cov}_{\mu_1}(X).
 \tag{7.1}
\]

Thus (DMMSE) is exactly

\[
 \mathbb E_B\,b^TA_1b\ge\kappa.                  \tag{7.2}
\]

The covariance SDE and the parallel derivative satisfy

\[
 dA_t=\sum_iH_{i,t}\,dB_{i,t}-A_t^2dt,
 \qquad
 {d\over dt}M_t=A_tM_t,
 \qquad M_0=I.                                     \tag{7.3}
\]

Therefore

\[
 \boxed{
 J_t:=A_tM_t,qquad
 dJ_t=\sum_iH_{i,t}M_t\,dB_{i,t},qquad
 \mathbb EJ_t=I.}
 \tag{7.4}
\]

The martingale \(J_t\) is the transpose version of the endpoint pairing
appearing in \(M_t^Tc_t\).

### 7.1 The MMSE floor closes the near-linear branch directly

Let \(f\) be a normalized first eigenfunction with energy \(\lambda\), and
write

\[
 f=a\cdot x+r,qquad
 \|r\|_{L^2(\mu)}=\delta,qquad
 |a|=\sqrt{1-\delta^2},qquad b={a\over|a|}.        \tag{7.5}
\]

At time one,

\[
 c_1:=\operatorname {Cov}_{\mu_1}(X,f)
 =A_1a+e_1,qquad
 e_1:=\operatorname {Cov}_{\mu_1}(X,r).            \tag{7.6}
\]

For each posterior and each unit \(u\),

\[
 |u\cdot e_1|^2
 \le\operatorname {Var}_1(u\cdot X)
      \operatorname {Var}_1(r)
 \le\operatorname {Var}_1(r),                     \tag{7.7}
\]

because \(A_1\preceq I\).  Taking the supremum in \(u\) and averaging
gives

\[
 \boxed{\mathbb E|e_1|^2\le\delta^2.}             \tag{7.8}
\]

Similarly,

\[
 |c_1|^2\le\operatorname {Var}_1(f)
 \le\mathbb E_1|\nabla f|^2,                      \tag{7.9}
\]

where the last step is the posterior Poincare inequality
\(C_P(\mu_1)\le1\).  Averaging the tilts back to \(\mu\) yields

\[
 \boxed{\mathbb E|c_1|^2\le\lambda.}              \tag{7.10}
\]

The MMSE floor has a stronger \(L^2\) consequence.  Pointwise,

\[
 |A_1b|\ge b^TA_1b\ge0.
\]

Hence Jensen and (7.2) give

\[
 \|A_1b\|_{L^2}
 \ge\mathbb E|A_1b|
 \ge\mathbb E\,b^TA_1b
 \ge\kappa.                                       \tag{7.11}
\]

The reverse triangle inequality in \(L^2\), followed by
(7.8)--(7.11), now gives

\[
\begin{aligned}
 \sqrt\lambda
 &\ge\|c_1\|_{L^2}\\
 &=\|\,|a|A_1b+e_1\|_{L^2}\\
 &\ge |a|\|A_1b\|_{L^2}-\|e_1\|_{L^2}\\
 &\ge \kappa\sqrt{1-\delta^2}-\delta.             \tag{7.12}
\end{aligned}
\]

This is (0.4).  In particular, if \(0<\kappa\le1\) and
\(\delta\le\kappa/4\), then the right side is at least \(\kappa/2\), so

\[
 \boxed{\lambda\ge\kappa^2/4.}                    \tag{7.13}
\]

No use of \(M_t\), \(J_t\), or final alignment is needed.  The endpoint
covariance relation alone turns (DMMSE) into the desired near-linear
spectral lower bound.

### 7.2 Why it still does not imply final alignment

In the exactly affine case \(e_1=0\), the proposed alignment integrand in
the selected direction is

\[
 {\langle M_1b,A_1b\rangle^2\over|A_1b|^2}
 ={(b^TJ_1b)^2\over b^TA_1^2b}.                    \tag{7.14}
\]

The floor gives

\[
 \mathbb E\,b^TA_1b\ge\kappa,
 \qquad
 \mathbb E\,b^TA_1^2b\ge\kappa^2,                \tag{7.15}
\]

where the second inequality follows from
\(b^TA_1^2b\ge(b^TA_1b)^2\) and Jensen.  But an average lower bound on the
denominator does not control the expectation of its reciprocal weighted by
\((b^TJ_1b)^2\).  The martingale identity supplies
\(\mathbb E(b^TJ_1b)=1\), which points in the wrong direction for an upper
bound on (7.14).

Here is a one-time algebraic model showing that all the natural terminal
first moments and trace scales are insufficient.  Fix \(b=e_1\), let
\(p=1/n\), and on the selected coordinate set

\[
\begin{array}{c|ccc}
 &M_1&J_1&A_1=J_1/M_1\\ \hline
 \text{rare event }(p)&n&1/2&1/(2n)\\[2mm]
 \text{common event }(1-p)&n/(n-1)&(2n-1)/(2(n-1))&(2n-1)/(2n).
\end{array}                                        \tag{7.16}
\]

On every orthogonal coordinate put

\[
 M_1=2,qquad J_1=1,qquad A_1={1\over2}.           \tag{7.17}
\]

All matrices are diagonal, positive, and satisfy \(J_1=A_1M_1\).  Direct
calculation gives

\[
 \mathbb EJ_1=I,qquad
 \mathbb EM_1=2I,qquad
 \mathbb EA_1\succeq {1\over2}I.                  \tag{7.18}
\]

Moreover,

\[
 \mathbb E\operatorname {Tr}(M_1^TM_1)
 =4(n-1)+n+{n\over n-1}=O(n).                      \tag{7.19}
\]

Nevertheless, if \(a=b\) and \(c_1=A_1a\), then the quotient (7.14) is
\(M_{1,11}^2\), and

\[
 \mathbb E{\langle M_1b,c_1\rangle^2\over|c_1|^2}
 =n+{n\over n-1}\asymp n.                         \tag{7.20}
\]

This model is not claimed to solve the localization SDE.  It shows exactly
what the terminal data fail to rule out: a rare event with small posterior
variance and large parallel amplification.  Any proof of final alignment
from the actual flow must use additional pathwise or quadratic-variation
structure that excludes this correlation.

With a residual \(e_1\ne0\), the quotient is even less stable on events
where \(|A_1b|\) is small, because a globally small \(L^2\) error can rotate
or nearly cancel \(A_1a\) on rare events.  The direct norm argument
(7.12) is robust to this; the reciprocal-weighted alignment quotient is
not.

---

## 8. KLS ledger and the precise live target

The logical status of every potentially dangerous step is:

1. **Unconditional channel calculus:**

   \[
   Dm=A,\quad 0\preceq A\preceq I,\quad
   \mathbb EA=I-J(\mu*\gamma).
   \]

2. **Unconditional entropy conclusion at the input entropy scale:**

   \[
   \operatorname {Tr}\mathbb EA
   \ge {n\over2}N(X)
   \ge {n\over4\pi eL_\infty(X)^2}.
   \]

3. **[SLICING INPUT]:**

   Replacing \(N(X)\) by a universal constant is the
   slicing/hyperplane conjecture.  KLS implies slicing, but this is weaker
   than the direct KLS input below and is not by itself directional.

4. **[KLS INPUT]:**

   \[
   C_P(\mu)\le C
   \quad\Longrightarrow\quad
   \mathbb EA\succeq {1\over C+2}I.
   \]

5. **[KLS CONSEQUENCE, NOT KNOWN EQUIVALENCE]:** (DMMSE) is implied by KLS,
   and failure of (DMMSE) would refute KLS.  No reduction from (DMMSE) to
   the full Poincare inequality was found.

6. **[OPEN OPERATOR STEP]:** Prove directly, for the special heat-smoothed
   density \(q=\mu*\gamma\),

   \[
   J(q)\preceq(1-\kappa)I,                         \tag{8.1}
   \]

   using its log-concavity, covariance \(2I\), and curvature bound

   \[
   0\preceq-D^2\log q=I-A\preceq I.                \tag{8.2}
   \]

   Equation (8.1) is not a consequence of the trace entropy estimate and
   should not be inserted as a generic Fisher inequality: it is exactly
   the desired directional theorem.

7. **Near-linear payoff:** If the open operator step is proved, (7.12)
   immediately rules out a near-linear small-gap eigenfunction.  Final
   alignment is not required for this implication.

The narrowest useful theorem exposed by the audit is therefore (DMMSE)
itself, or even its version only for the data-dependent direction
\(b=a/|a|\) selected by a near-linear first eigenfunction.  The latter
restriction might permit use of the eigenfunction/Bochner structure and is
strictly more targeted than proving a full operator bound for every
direction and every isotropic log-concave input.

## Reference used for the entropy step

* M. Marsiglietti and V. Kostina, [*A lower bound on differential entropy
  of log-concave random vectors with
  applications*](https://arxiv.org/abs/1704.07766).  The vector constants
  are dimension-dependent in general and become dimension-free under
  additional hypotheses; a dimension-free general bound would enter the
  slicing problem.
