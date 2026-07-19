# Upper-curvature one-form coercivity: exact equivalence and surviving branch

## 0. Verdict

Let

\[
 d\mu=Z^{-1}e^{-V}dx,
 \qquad 0\preceq H:=D^2V\preceq2I,
 \qquad \operatorname {Cov}(\mu)=I.
\]

The proposed estimate

\[
 \int|w|^2d\mu
 \le C\left\{\int\|Dw\|_{\rm HS}^2d\mu
                  +\int w^THw\,d\mu\right\},
 \qquad w=\nabla f,                                  \tag{OF}
\]

would prove the fixed-Gaussian subclass and hence, by the posterior reverse
theorem, KLS.  It is not merely sufficient: its optimal constant is exactly
the scalar Poincare constant of \(\mu\).  Thus (OF), even after restricting
to gradient fields and even under the displayed upper-curvature hypotheses,
is the remaining KLS-equivalent assertion rather than a consequence of
Bochner.

There are two genuine dimension-free branches:

1. the average-curvature identity gives (OF) for gradients whose variation
   around a constant vector is a sufficiently small fraction of their total
   mass; and
2. upper-curvature score rigidity gives a numerical gap for bottom spectral
   vectors sufficiently close in \(L^2\) to the affine space.

A small bottom spectral value automatically avoids the first branch, and
the existing score argument forces it to avoid the second.  The survivor is
an almost-centered, genuinely nonlinear exact field which varies slowly and
aligns with the low-eigenvalue space of \(H(x)\) adaptively.  Controlling
that survivor is exactly the unresolved nonlinear branch.

## 1. Exact scalar/one-form equivalence

Let

\[
 L=\Delta-\nabla V\cdot\nabla,
 \qquad A=-L
\]

be the nonnegative self-adjoint generator on \(L^2(\mu)\), and let
\(\lambda=\lambda _1(\mu)\).  For a smooth core function \(f\), the
integrated Bochner identity is

\[
 \|Af\|_2^2
 =\int\|D^2f\|_{\rm HS}^2d\mu
  +\int\nabla f^TH\nabla f\,d\mu.                    \tag{1.1}
\]

Taking \(w=\nabla f\), the right side of (OF) is therefore exactly
\(\|Af\|_2^2\), while

\[
 \int|w|^2d\mu=\langle f,Af\rangle.                 \tag{1.2}
\]

Constants may be subtracted from \(f\) without changing either quantity.
The spectral theorem on \(\mathbf1^\perp\) gives

\[
 \sup_{0\ne f\in\operatorname {Dom}(A)\cap\mathbf1^\perp}
 \frac{\langle f,Af\rangle}{\|Af\|_2^2}
 =\frac1\lambda=C_P(\mu).                            \tag{1.3}
\]

Indeed, the upper bound follows from \(A\succeq\lambda I\).  For the
reverse bound, take a unit vector in the spectral subspace
\(\mathbf1_{[\lambda,\lambda+\varepsilon]}(A)\) and send
\(\varepsilon\downarrow0\).  Thus (1.3) also covers a non-attained or
continuous spectral edge.  Core closure extends (1.1)--(1.3) to the exact
one-form form domain.  No upper bound on \(H\) was used in this equivalence.

Equivalently, on curl-free fields the Witten one-form quadratic form is the
square of the weighted divergence:

\[
 \int\|Dw\|_{\rm HS}^2+\int w^THw
 =\int|\operatorname {div}w-\nabla V\cdot w|^2d\mu,
 \qquad w=\nabla f.                                  \tag{1.4}
\]

Consequently calling (OF) a one-form, Korn, or divergence coercivity
estimate does not weaken its spectral content.

## 2. What average upper curvature proves

Integration by parts and the matrix Cramer--Rao inequality give

\[
 \mathbb EH=\mathbb E[\nabla V\nabla V^T]\succeq I. \tag{2.1}
\]

For completeness,
\(\mathbb E[X\nabla V^T]=I\), and positivity of the covariance matrix of
\((X,\nabla V(X))\) gives (2.1).  The upper Hessian hypothesis also gives
\(\mathbb EH\preceq2I\).

Let \(w=\nabla f\), put \(m=\mathbb Ew\), and write \(u=w-m\).  From
(2.1), \(H\preceq2I\), and the triangle inequality in the
\(H(x)\)-seminorm,

\[
 \begin{aligned}
 |m|^2
 &\le\mathbb E[m^THm]\\
 &\le2\mathbb E[w^THw]+2\mathbb E[u^THu]\\
 &\le2\mathbb E[w^THw]+4\mathbb E|u|^2.
 \end{aligned}                                       \tag{2.2}
\]

Hence the following is a genuine global, dimension-free branch:

\[
 \boxed{\quad
 \mathbb E|w|^2
 \le2\mathbb E[w^THw]+5\mathbb E|w-\mathbb Ew|^2.
 \quad}                                               \tag{2.3}
\]

In particular, if

\[
 \mathbb E|w-\mathbb Ew|^2\le\eta\mathbb E|w|^2,
 \qquad 0\le\eta<\frac15,
\]

then

\[
 \boxed{\quad
 \mathbb E|w|^2
 \le\frac2{1-5\eta}\mathbb E[w^THw].
 \quad}                                               \tag{2.4}
\]

Thus average curvature really does control a near-constant gradient; no
Poincare inequality occurs in (2.2)--(2.4).

The tempting next step is to control the centered field \(u\) by \(Du\).
Define the optimal exact-field constant

\[
 K_\mu=\sup_{f}
 \frac{\mathbb E|\nabla f-\mathbb E\nabla f|^2}
      {\mathbb E\|D^2f\|_{\rm HS}^2}.                \tag{2.5}
\]

Componentwise Poincare gives \(K_\mu\le C_P(\mu)\).  Conversely, a unit
spectral-window vector in
\([\lambda,\lambda+\varepsilon]\) satisfies

\[
 \mathbb E|\nabla f|^2\ge\lambda,
 \qquad
 |\mathbb E\nabla f|\le\|Af\|_2\le\lambda+\varepsilon,
 \qquad
 \mathbb E\|D^2f\|_{\rm HS}^2\le(\lambda+\varepsilon)^2.
\]

Here the middle inequality follows by testing \(Af\) against the isotropic
coordinate functions, and the last follows from (1.1).  Sending
\(\varepsilon\downarrow0\) gives

\[
 \boxed{\qquad C_P(\mu)-1\le K_\mu\le C_P(\mu).\qquad} \tag{2.6}
\]

Therefore inserting a universal centered-gradient/Korn estimate into
(2.3) is circular up to the exact additive constant one.

## 3. Bottom spectral windows and the precise survivor

Let \(f_\varepsilon\) be a centered unit vector in the spectral window
\([\lambda,\lambda+\varepsilon]\), and set

\[
 q_\varepsilon=\mathbb E|\nabla f_\varepsilon|^2,
 \qquad
 r_\varepsilon=\|Af_\varepsilon\|_2^2,
 \qquad
 W_\varepsilon=\frac{\nabla f_\varepsilon}
                      {\sqrt{q_\varepsilon}}.
\]

Then

\[
 \mathbb E|W_\varepsilon|^2=1,
 \qquad
 \mathbb E\|DW_\varepsilon\|_{\rm HS}^2
 +\mathbb E[W_\varepsilon^THW_\varepsilon]
 =\frac{r_\varepsilon}{q_\varepsilon}
 \le\frac{(\lambda+\varepsilon)^2}{\lambda},         \tag{3.1}
\]

and isotropy gives

\[
 |\mathbb EW_\varepsilon|^2
 \le\frac{r_\varepsilon}{q_\varepsilon}
 \le\frac{(\lambda+\varepsilon)^2}{\lambda}.         \tag{3.2}
\]

Thus, as \(\varepsilon=o(\lambda)\), a hypothetical small \(\lambda\)
produces a unit exact field which is almost centered and has both deformation
and curvature energy \(O(\lambda)\).  It lies at distance
\(1-O(\lambda)\) from the near-constant branch (2.4).

There is nevertheless a separate, proved function-specific theorem.  If
\(\delta_\varepsilon\) denotes the \(L^2\)-distance of
\(f_\varepsilon\) from the affine space, upper-curvature score rigidity
gives, for \(H\preceq\beta I\),

\[
 (1-\lambda)\sqrt{1-\delta_\varepsilon^2}
 \le\sqrt{\beta-1}\,\delta_\varepsilon+\varepsilon.  \tag{3.3}
\]

For the fixed-Gaussian value \(\beta=2\), a sequence with
\(\delta_\varepsilon\le1/2\) forces

\[
 \lambda\ge1-\frac1{\sqrt3}.                         \tag{3.4}
\]

Equations (2.4) and (3.3) are genuine global coercivity statements.  What
they do not control is a bottom window with an almost-centered gradient and
a scalar test function bounded away from the affine space.  Applying (OF)
to that survivor would simply assume the desired conclusion through (1.3).

## 4. Extra structure of an actual fixed-Gaussian output

For

\[
 S=\frac{X+G}{\sqrt2}
\]

let \(A_y=\operatorname {Cov}(X\mid X+G=y)\).  The output potential obeys

\[
 H(s)=2(I-A_{\sqrt2s}).                               \tag{4.1}
\]

Consequently the curvature part of a low-energy field is

\[
 \mathbb E[W^THW]
 =2\mathbb E\left(|W(S)|^2
       -W(S)^TA_{\sqrt2S}W(S)\right).                 \tag{4.2}
\]

Since \(0\preceq A_y\preceq I\), a survivor from (3.1) must choose, at
most output points, a direction in which the strongly log-concave posterior
nearly saturates its covariance bound.  Simultaneously \(DW\) is small, so
those locally saturating directions must form a slowly varying curl-free
field.

This is more structure than the abstract assumptions
\(0\preceq H\preceq2I\), but the current weighted-Fisher estimate controls
posterior error only in each fixed deterministic direction.  It does not
control an adaptive direction \(W(S)\).  A valid new route would therefore
be an adaptive posterior-saturation coherence theorem: near saturation in
(4.2), together with small \(DW\), must force a coherent Gaussian factor;
the almost-centered condition (3.2) would then exclude that factor.  No such
theorem is proved here, and restating it as (OF) would be circular.

## 5. Stress tests

### 5.1 Isotropic interval and Gamma products

Fixed Gaussian regularization acts coordinatewise, so products stay
products and their gaps tensorize.  For the isotropic interval, the input
gap is

\[
 a_I=\frac{\pi^2}{12}.
\]

The proved forward two-law inequality gives for its normalized
unit-Gaussian output

\[
 b_I\ge\frac{\sqrt{9+16a_I}-3}{2}>0.85.              \tag{5.1}
\]

For the standardized Gamma law of shape \(k\ge1\), the exact input gap is

\[
 a_k=\frac{k^2}{(k+1)^2}.
\]

Thus every fixed-Gaussian output satisfies

\[
 b_k\ge\frac{\sqrt{9+16a_k}-3}{2}
 \ge\frac{\sqrt{13}-3}{2}>0.30.                      \tag{5.2}
\]

For arbitrary tensor powers the optimal constant in (OF) is exactly
\(1/\min_i b_i\), by (1.3), so it stays below \(3.31\).  Slow
one-coordinate tests attain the relevant tensor edge.  Hence interval,
one-sided exponential, and Gamma products exhibit neither a hidden factor
of the dimension nor any gain beyond ordinary tensorized spectral
coercivity.

### 5.2 Rotating nullspaces and radial models

The matrix average (2.1) alone cannot prevent an adiabatically rotating
kernel.  On an abstract interval of length \(L\), let

\[
 H(t)=2v(t)v(t)^T,
 \quad v(t)=(\cos(\pi t/L),\sin(\pi t/L)),
 \quad w(t)=(-\sin(\pi t/L),\cos(\pi t/L)).
\]

Then the uniform average of \(H\) is \(I\), while
\(w^THw=0\) and \(|w'|^2=\pi^2/L^2\).  This is not a counterexample to
(OF): the matrix field is not supplied as the Hessian of the base density,
and the rotating vector field is not a scalar gradient in the corresponding
two-dimensional space.  The toy model identifies exactly why both Hessian
integrability and curl-freeness must be used by any proof.

The natural curl-free radial realization pays the missing angular energy.
For any isotropic radial law in \(\mathbb R^n\), take \(f(x)=|x|\).  Away
from the origin,

\[
 w=e_r,
 \qquad Dw=\frac{I-e_re_r^T}{|x|},
 \qquad \|Dw\|_{\rm HS}^2=\frac{n-1}{|x|^2}.
\]

Since \(\mathbb E|X|^2=n\), Jensen gives

\[
 \mathbb E\|Dw\|_{\rm HS}^2
 =(n-1)\mathbb E|X|^{-2}
 \ge\frac{n-1}{n}.                                   \tag{5.3}
\]

Thus even if the radial direction is a pointwise curvature null direction,
its rotation costs order one.  At the opposite extreme a constant field
pays at least its full norm through (2.1).  A bad field must interpolate
between these extremes in a genuinely nonradial and adaptive fashion.

### 5.3 The exponential cone

The cone model in `gaussian_curvature_korn.md` has conditional slices
\(q_a={\rm Unif}[-a,a]*\gamma\).  Its exact conditional Poisson field has
deformation energy \(O(a^{-2})\) and boundary-layer curvature energy of
order at least \(a^{-1}\), while its field norm stays of order one for
large \(a\).  The last assertion follows from

\[
 \mathbb E|F_a|^2
 =a^{-2}\operatorname {Var}(\mathbb E[Y\mid Y+G])
 \ge\frac13-\frac1{a^2},
\]

because \(\operatorname {Var}(Y)=a^2/3\) and the posterior variance is at
most one.  These estimates show that the curvature term cannot be discarded
or bounded by deformation energy slice by slice.  They do **not** refute
(OF), where that curvature term occurs with the favorable sign on the
right-hand side.

It does not produce a global survivor.  The scalar width has the
standardized Gamma\((2,1)\) marginal, whose exact gap is \(4/9\); localizing
a joint gradient in the far-width tail therefore incurs a universal scalar
derivative cost.  The slices with large local ratio also have exponentially
small marginal weight.  Full Gaussian regularization puts the joint law in
the class \(H\preceq2I\), but does not turn the failed slice-local estimate
into a proof of (OF).  The cone stress test therefore confirms that any
successful argument must be global in the scalar and transverse variables;
it neither refutes nor proves the KLS-equivalent global estimate.

## 6. Audit conclusion

The upper-curvature work now supplies three valid pieces:

1. fixed-Gaussian outputs are a quantitatively equivalent KLS target by
   posterior reverse smoothing;
2. the pointwise WFI estimate gives a dimension-free directional MMSE floor;
3. average curvature and score rigidity close near-constant-gradient and
   near-affine bottom branches, respectively.

The proposed global exact-field coercivity (OF) is not a fourth piece: by
(1.3) it is exactly the missing spectral gap.  Products, radial rotation,
and the cone do not contradict it, but they show that fixed directions,
pointwise average curvature, and slice-local propagation are insufficient.
The materially narrower unresolved mechanism is adaptive posterior
saturation coherence in (4.2), not another generic one-form Poincare
estimate.
