# Convolution stability inverse: the affine tangent obstruction and a corrected defect

## 1. Verdict

Let \(X,Y\) be iid with centered isotropic log-concave law \(\mu\), let

\[
 S=\frac{X+Y}{\sqrt2},\qquad \nu=\mathcal L(S),
\]

and define

\[
 Kf(x)=\mathbb E_Y f\!\left(\frac{x+Y}{\sqrt2}\right).
\]

The proposed estimates

\[
 2\operatorname{Var}_\mu(Kf)
 \le \int |\nabla f|^2\,d\nu+C\|r_f\|_2^2
 \tag{1.1}
\]

and

\[
 \mathbb E(\Delta_4f)^2
 \ge c\left(
 \operatorname{Var}_\nu f-\int|\nabla f|^2\,d\nu
 \right)_+
 \tag{1.2}
\]

are false for every non-Gaussian \(\mu\), even in dimension one and even for a smooth, strongly log-concave, exactly isotropic density.

The obstruction is first order.  The Hoeffding residual vanishes on affine functions and changes by \(O(t)\) under a perturbation \(\ell+tg\), so its squared norm is \(O(t^2)\).  Unless \(\nu\) is Gaussian, the excess on the right of (1.2), and the left-minus-first-term in (1.1), can change by \(O(t)\).  Choosing the sign of \(t\) disproves any finite constant.

More precisely:

* Either estimate, if valid for all \(f\), forces the Gaussian Stein identity for \(\nu\), hence forces both \(\nu\) and \(\mu\) to be Gaussian.
* An explicit smooth strongly log-concave counterexample is given below by a quadratic tilt of a Gaussian and the cubic perturbation \(f_\varepsilon(s)=s-\varepsilon s^3\).  The ratio in (1.1) diverges like \(8/(129\varepsilon)\), while the inverse ratio in (1.2) diverges like \(2/(129\varepsilon)\).
* Uniform regular simplices fail already with a linear-plus-quadratic test because their vertex-direction third moment is nonzero.  Smooth strictly convex approximations inherit the failure.
* On pure quadratic modes, either estimate would imply the generalized quadratic variance bound
  \[
  \operatorname{Var}(X^TAX)\lesssim\operatorname{tr}(A^2).
  \]
  Taking \(A=I\) gives the thin-shell variance conjecture.  Thus removing the affine tangent counterexample still leaves a KLS-strength open obstruction.

The minimal repair is to add an affine--Stein cross defect.  If

\[
 a_f=\mathbb E_\nu[S(f-\mathbb E f)],
 \qquad
 \tau_f=a_f-\mathbb E_\nu\nabla f,
 \tag{1.3}
\]

then the missing first-order term is exactly \(2\langle a_f,\tau_f\rangle\).  Subtracting the \(L^2(\nu)\)-affine projection makes this identity exact.  A corrected stability statement must either include this term or restrict to the affine-orthogonal sector; the latter still contains the generalized quadratic variance problem.

## 2. Canonical centering and equivalence of the two proposed forms

Put

\[
 F(x,y)=f\!\left(\frac{x+y}{\sqrt2}\right),
 \qquad m=\mathbb E_\nu f,
\]

and use the canonical Hoeffding decomposition

\[
 h_f(x)=Kf(x)-m,
 \qquad
 r_f(x,y)=F(x,y)-m-h_f(x)-h_f(y).
 \tag{2.1}
\]

Equivalently,

\[
 r_f(x,y)=F(x,y)+m-Kf(x)-Kf(y).
\]

If \(f\) is centered, this is simply \(F-Kf(x)-Kf(y)\), so it agrees with the residual in the question after harmless centering.  One has

\[
 \operatorname{Var}_\nu f
 =2\|h_f\|_{L^2(\mu)}^2+\|r_f\|_{L^2(\mu^2)}^2
 \tag{2.2}
\]

and

\[
 \mathbb E(\Delta_4f)^2=4\|r_f\|_2^2.
 \tag{2.3}
\]

Define

\[
 \mathcal D(f)
 :=\operatorname{Var}_\nu f-\int|\nabla f|^2\,d\nu,
 \qquad
 \mathcal Q(f)
 :=2\operatorname{Var}_\mu(Kf)-\int|\nabla f|^2\,d\nu.
 \tag{2.4}
\]

Then (2.2) gives the exact relation

\[
 \boxed{\mathcal D(f)=\mathcal Q(f)+\|r_f\|_2^2.}
 \tag{2.5}
\]

Consequently, (1.1) with constant \(C\ge0\) implies (1.2) with

\[
 c=\frac4{C+1}.
\]

Conversely, (1.2) implies (1.1), for example with \(C=4/c\).  Thus the two statements have the same obstruction.

## 3. General no-go theorem: the estimate characterizes the Gaussian

For \(a\in\mathbb R^n\), let

\[
 \ell_a(z)=a\cdot z.
\]

Its Hoeffding pieces are

\[
 h_{\ell_a}(x)=\frac{a\cdot x}{\sqrt2},
 \qquad r_{\ell_a}=0,
\]

and isotropicity gives

\[
 \mathcal D(\ell_a)=\mathcal Q(\ell_a)=0.
\]

For a smooth \(g\), define the Stein bilinear form

\[
 \mathcal B_\nu(a,g)
 :=\mathbb E_\nu[(a\cdot S)g(S)]
   -\mathbb E_\nu[a\cdot\nabla g(S)].
 \tag{3.1}
\]

### Lemma 3.1 (first variation at the affine space)

For every \(a,g\) and scalar \(t\),

\[
 \boxed{
 \begin{aligned}
 \mathcal D(\ell_a+tg)
 &=2t\,\mathcal B_\nu(a,g)+t^2\mathcal D(g),\\
 \mathcal Q(\ell_a+tg)
 &=2t\,\mathcal B_\nu(a,g)+t^2\mathcal Q(g),\\
 r_{\ell_a+tg}&=t\,r_g.
 \end{aligned}}
 \tag{3.2}
\]

#### Proof

The first identity follows by polarizing variance and Dirichlet energy.  For the second, symmetry in \(X,Y\) gives

\[
 \mathbb E_\nu[(a\cdot S)g(S)]
 =\sqrt2\,\mathbb E[(a\cdot X)g(S)].
\]

Since \(h_{\ell_a}(X)=a\cdot X/\sqrt2\),

\[
 2\langle h_{\ell_a},h_g\rangle_{L^2(\mu)}
 =\mathbb E_\nu[(a\cdot S)g(S)].
\]

Polarizing \(\mathcal Q\) now gives the same cross term (3.1).  The last identity follows because the Hoeffding residual is linear and annihilates affine functions.  \(\square\)

### Theorem 3.2 (Gaussian characterization)

Suppose either (1.1) with a finite \(C\), or (1.2) with a positive \(c\), holds for every smooth \(f\).  Then

\[
 \nu=N(0,I)
\qquad\text{and}\qquad
 \mu=N(0,I).
 \tag{3.3}
\]

Conversely, for Gaussian \(\mu\), (1.1) holds with \(C=0\).

#### Proof

If \(\mathcal B_\nu(a,g)\ne0\), choose the sign of \(t\) so that

\[
 t\mathcal B_\nu(a,g)>0.
\]

By (3.2), either \(\mathcal D(\ell_a+tg)\) or \(\mathcal Q(\ell_a+tg)\), as appropriate, is positive of order \(|t|\), whereas \(\|r_{\ell_a+tg}\|_2^2\) and \(\mathbb E(\Delta_4(\ell_a+tg))^2\) are of order \(t^2\).  This contradicts either proposed inequality as \(t\to0\).  Therefore

\[
 \mathbb E[S_i g(S)]=\mathbb E[\partial_i g(S)]
 \tag{3.4}
\]

for every \(i\) and every smooth \(g\).  Applying (3.4) to Fourier exponentials, by approximation, gives

\[
 \partial_{t_i}\widehat\nu(t)=-t_i\widehat\nu(t).
\]

Since \(\widehat\nu(0)=1\), this yields

\[
 \widehat\nu(t)=e^{-|t|^2/2}.
\]

Also

\[
 \widehat\nu(t)=\widehat\mu(t/\sqrt2)^2.
\]

The right side never vanishes.  Its continuous square root equals \(1\) at the origin, so

\[
 \widehat\mu(u)=e^{-|u|^2/2}.
\]

Thus both laws are Gaussian.

Conversely, if \(\mu\) is Gaussian, so is \(\nu\), and Gaussian Poincare together with (2.2) gives

\[
 2\operatorname{Var}_\mu(Kf)
 =\operatorname{Var}_\nu f-\|r_f\|_2^2
 \le\int|\nabla f|^2\,d\nu-\|r_f\|_2^2.
\]

This is stronger than (1.1) with \(C=0\).  \(\square\)

This theorem is distributional: log-concavity is not used.  Log-concavity therefore cannot repair the missing first variation.

## 4. An explicit smooth strongly log-concave counterexample

The preceding theorem already proves failure, but the following example gives exact constants.

Let \(\phi\) be the standard Gaussian density and, for \(0<\alpha<1/2\), let \(W_\alpha\) have density

\[
 p_\alpha(w)
 =\frac{1+\alpha w^2}{1+\alpha}\,\phi(w).
 \tag{4.1}
\]

Its log-density has second derivative

\[
 -1+\frac{2\alpha(1-\alpha w^2)}
 {(1+\alpha w^2)^2}
 \le -1+2\alpha<0.
 \tag{4.2}
\]

Hence it is smooth and strongly log-concave.  Set

\[
 \sigma_\alpha^2
 =\mathbb EW_\alpha^2
 =\frac{1+3\alpha}{1+\alpha},
 \qquad
 X_\alpha=\frac{W_\alpha}{\sigma_\alpha}.
 \tag{4.3}
\]

Then \(\mu_\alpha=\mathcal L(X_\alpha)\) is centered and exactly isotropic.  Its potential, up to an additive constant, is

\[
 V_\alpha(x)
 =\frac{\sigma_\alpha^2x^2}{2}
  -\log(1+\alpha\sigma_\alpha^2x^2),
 \tag{4.4}
\]

and (4.2) shows \(V_\alpha''>0\) uniformly.

The standardized fourth and sixth moments are

\[
 m_4
 =\frac{(3+15\alpha)(1+\alpha)}
 {(1+3\alpha)^2},
 \qquad
 m_6
 =\frac{(15+105\alpha)(1+\alpha)^2}
 {(1+3\alpha)^3}.
 \tag{4.5}
\]

In particular,

\[
 m_4-3
 =-\frac{12\alpha^2}{(1+3\alpha)^2}<0.
 \tag{4.6}
\]

For any centered symmetric variance-one law, put

\[
 s_4=\mathbb ES^4=\frac{m_4+3}{2},
 \qquad
 s_6=\mathbb ES^6=\frac{m_6+15m_4}{4}.
 \tag{4.7}
\]

For

\[
 f_\varepsilon(s)=s-\varepsilon s^3,
 \tag{4.8}
\]

one computes

\[
 Kf_\varepsilon(x)
 =\frac{x}{\sqrt2}
  -\frac{\varepsilon}{2\sqrt2}(x^3+3x),
 \tag{4.9}
\]

and

\[
 r_{f_\varepsilon}(x,y)
 =-\frac{3\varepsilon}{2\sqrt2}
 \left[y(x^2-1)+x(y^2-1)\right].
 \tag{4.10}
\]

Symmetry and independence give

\[
 \|r_{f_\varepsilon}\|_2^2
 =\frac94(m_4-1)\varepsilon^2,
 \qquad
 \mathbb E(\Delta_4f_\varepsilon)^2
 =9(m_4-1)\varepsilon^2.
 \tag{4.11}
\]

Direct moment expansion gives

\[
 \mathcal D(f_\varepsilon)
 =(3-m_4)\varepsilon
 +(s_6-9s_4)\varepsilon^2.
 \tag{4.12}
\]

Thus the positive excess is linear in \(\varepsilon\), while both residual defects are quadratic.

For a completely rational instance, take \(\alpha=1/4\).  Then

\[
 \sigma_\alpha^2=\frac75,\qquad
 m_4=\frac{135}{49},\qquad
 m_6=\frac{4125}{343},
 \tag{4.13}
\]

The corresponding potential is, up to an additive constant,

\[
 V(x)=\frac7{10}x^2-\log\left(1+\frac7{20}x^2\right),
 \qquad V''(x)\ge\frac7{10}.
 \tag{4.14}
\]

The exact defect formulas become

\[
 \begin{aligned}
 \mathcal D(f_\varepsilon)
 &=\frac{12}{49}\varepsilon
   -\frac{4308}{343}\varepsilon^2,\\
 \|r_{f_\varepsilon}\|_2^2
 &=\frac{387}{98}\varepsilon^2,\\
 \mathcal Q(f_\varepsilon)
 &=\frac{12}{49}\varepsilon
   -\frac{11325}{686}\varepsilon^2,\\
 \mathbb E(\Delta_4f_\varepsilon)^2
 &=\frac{774}{49}\varepsilon^2.
 \end{aligned}
 \tag{4.15}
\]

Consequently,

\[
 \frac{\mathcal Q(f_\varepsilon)}
 {\|r_{f_\varepsilon}\|_2^2}
 =\frac{8}{129}\frac1\varepsilon+O(1)
 \longrightarrow\infty,
 \tag{4.16}
\]

and

\[
 \frac{\mathcal D(f_\varepsilon)}
 {\mathbb E(\Delta_4f_\varepsilon)^2}
 =\frac{2}{129}\frac1\varepsilon+O(1)
 \longrightarrow\infty.
 \tag{4.17}
\]

This is an exactly isotropic, smooth, full-support, uniformly convex counterexample in dimension one.  Products give the same counterexample in every dimension by using a function of one coordinate.

## 5. Simplex and skew-mode counterexamples

There is an even lower-degree obstruction whenever a one-dimensional marginal has nonzero third moment.

Let \(u\) be a unit vector, put

\[
 Z=u\cdot X,\qquad m_3=\mathbb EZ^3,\qquad m_4=\mathbb EZ^4,
\]

and consider

\[
 f_t(s)=u\cdot s+t\big((u\cdot s)^2-1\big).
 \tag{5.1}
\]

Writing \(z=u\cdot x\), one has exactly

\[
 Kf_t(x)=\frac z{\sqrt2}+\frac t2(z^2-1),
 \qquad
 r_{f_t}(x,y)=t(u\cdot x)(u\cdot y).
 \tag{5.2}
\]

Since the third moment of \(u\cdot S\) is \(m_3/\sqrt2\),

\[
 \begin{aligned}
 \|r_{f_t}\|_2^2&=t^2,\\
 \mathbb E(\Delta_4f_t)^2&=4t^2,\\
 \mathcal Q(f_t)
 &=\sqrt2\,m_3t+
 \left(\frac{m_4-1}{2}-4\right)t^2,\\
 \mathcal D(f_t)
 &=\sqrt2\,m_3t+\frac{m_4-7}{2}t^2.
 \end{aligned}
 \tag{5.3}
\]

If \(m_3\ne0\), choosing the sign of \(t\) makes both ratios diverge like \(1/|t|\).

For the uniform law on an isotropic regular \(n\)-simplex, choose \(u\) pointing toward a vertex.  The corresponding barycentric coordinate has the \(\operatorname{Beta}(1,n)\) law.  After centering and variance normalization, its third moment is

\[
 m_3
 =\frac{2(n-1)\sqrt{n+2}}
 {(n+3)\sqrt n},
 \tag{5.4}
\]

which is nonzero for every \(n\ge2\) and tends to \(2\).  Thus (5.1)--(5.3) give explicit simplex counterexamples.

Smooth strictly convex approximations retain them quantitatively.  For example, if \(K\) is the isotropic simplex, use

\[
 U_{a,\eta,\delta}(x)
 =a(\rho_\delta*d_K^2)(x)+\eta|x|^2,
 \tag{5.5}
\]

then whiten the covariance.  Sending \(a\to\infty\) and \(\eta,\delta\to0\) makes all moments in (5.3) converge to the simplex moments.  The nonzero linear coefficient in (5.3) persists.  The same observation applies to any skew thin cone: a nonzero directional third moment is already enough, so no delicate high-dimensional cone asymptotics are needed.

As a one-dimensional simplex-like limit, \(X=E-1\), with \(E\sim\operatorname{Exp}(1)\), has

\[
 m_3=2,\qquad m_4=9.
\]

For (5.1),

\[
 \mathcal Q(f_t)=2\sqrt2\,t,
 \qquad
 \mathcal D(f_t)=2\sqrt2\,t+t^2,
 \qquad
 \|r_{f_t}\|_2^2=t^2.
 \tag{5.6}
\]

The facet/exponential geometry therefore realizes the obstruction with no quadratic remainder in \(\mathcal Q\).

## 6. Exact quadratic-mode audit

Let \(A=A^T\) and define

\[
 f_A(z)=z^TAz-\operatorname{tr}A.
 \tag{6.1}
\]

No symmetry assumption on \(\mu\) is needed.  Isotropicity gives

\[
 Kf_A(x)=\frac12\big(x^TAx-\operatorname{tr}A\big),
 \qquad
 r_{f_A}(x,y)=x^TAy.
 \tag{6.2}
\]

Therefore

\[
 \boxed{
 \begin{aligned}
 2\operatorname{Var}_\mu(Kf_A)
 &=\frac12\operatorname{Var}_\mu(X^TAX),\\
 \int|\nabla f_A|^2\,d\nu
 &=4\operatorname{tr}(A^2),\\
 \|r_{f_A}\|_2^2
 &=\operatorname{tr}(A^2),\\
 \mathbb E(\Delta_4f_A)^2
 &=4\operatorname{tr}(A^2),\\
 \mathcal D(f_A)
 &=\frac12\operatorname{Var}_\mu(X^TAX)
   -3\operatorname{tr}(A^2).
 \end{aligned}}
 \tag{6.3}
\]

Thus (1.1) on quadratic modes is exactly the bound

\[
 \operatorname{Var}_\mu(X^TAX)
 \le (8+2C)\operatorname{tr}(A^2).
 \tag{6.4}
\]

Similarly, (1.2) implies

\[
 \operatorname{Var}_\mu(X^TAX)
 \le \left(6+\frac8c\right)\operatorname{tr}(A^2).
 \tag{6.5}
\]

Conversely, a generalized quadratic variance bound controls the proposed inequality on this polynomial sector with the corresponding constants.

Taking \(A=I\) gives

\[
 \operatorname{Var}|X|^2\lesssim n,
 \tag{6.6}
\]

the thin-shell variance conjecture.  For centrally symmetric \(\mu\), the quadratic \(f_A\) is \(L^2(\nu)\)-orthogonal to every affine function.  Hence projecting away the affine counterexample does not eliminate (6.4): the repaired statement still contains the generalized thin-shell problem and cannot be proved here without new KLS-strength input.

## 7. The minimal modified defect

The failure can be separated exactly from the remaining nonlinear problem.

Assume \(f\) is centered and set

\[
 a_f=\mathbb E_\nu[Sf(S)],
 \qquad
 b_f=\mathbb E_\nu\nabla f(S),
 \qquad
 \tau_f=a_f-b_f.
 \tag{7.1}
\]

The vector \(\tau_f\) is the Stein discrepancy of the test function against affine directions.  Let

\[
 g_f(s)=f(s)-a_f\cdot s.
 \tag{7.2}
\]

Because \(\nu\) is isotropic,

\[
 \mathbb E_\nu[Sg_f(S)]=0,
\]

so \(g_f\) is \(L^2(\nu)\)-orthogonal to the affine space.  Also

\[
 r_{g_f}=r_f.
 \tag{7.3}
\]

### Lemma 7.1 (exact affine--Stein decomposition)

One has

\[
 \boxed{
 \begin{aligned}
 \mathcal D(f)
 &=\mathcal D(g_f)+2\langle a_f,\tau_f\rangle,\\
 \mathcal Q(f)
 &=\mathcal Q(g_f)+2\langle a_f,\tau_f\rangle.
 \end{aligned}}
 \tag{7.4}
\]

#### Proof

For \(\mathcal D\), variance is orthogonal under (7.2), while

\[
 \int|\nabla f|^2
 =|a_f|^2+2\langle a_f,\mathbb E\nabla g_f\rangle
  +\int|\nabla g_f|^2.
\]

Since

\[
 \mathbb E\nabla g_f=b_f-a_f=-\tau_f,
\]

the first identity follows.  The second follows by polarizing \(\mathcal Q\) and using Lemma 3.1 together with \(\mathbb E[Sg_f]=0\).  \(\square\)

Thus the minimal homogeneous correction forced by the counterexamples is

\[
 \boxed{
 \mathfrak R_{\mathrm{mod}}(f)
 :=\mathbb E(\Delta_4f)^2
   +2\big|\langle a_f,\tau_f\rangle\big|.
 }
 \tag{7.5}
\]

Indeed, if one could prove the affine-orthogonal estimate

\[
 \mathcal D(g)_+
 \le C\,\mathbb E(\Delta_4g)^2,
 \qquad \mathbb E[Sg(S)]=0,
 \tag{7.6}
\]

then (7.4) would give

\[
 \mathcal D(f)_+
 \le C\,\mathbb E(\Delta_4f)^2
   +2\big|\langle a_f,\tau_f\rangle\big|.
 \tag{7.7}
\]

The smooth cubic example has \(|a_f|\asymp1\), \(|\tau_f|\asymp\varepsilon\), and \(\|r_f\|_2^2\asymp\varepsilon^2\), so the added term has exactly the missing order.

This correction is also relevant for actual Rayleigh-critical functions.  If

\[
 R(f)=\frac{\operatorname{Var}_\nu f}
 {\int|\nabla f|^2\,d\nu}
\]

is stationary under every affine perturbation, first variation gives

\[
 a_f=R(f)b_f.
 \tag{7.8}
\]

Therefore

\[
 \langle a_f,\tau_f\rangle
 =\frac{R(f)-1}{R(f)}|a_f|^2\ge0
 \qquad\text{when }R(f)\ge1.
 \tag{7.9}
\]

The affine--Stein term is consequently a genuine part of the excess of a top mode, not an artifact that criticality removes.

Estimate (7.6) is the honest residual-stability question after removing the fatal tangent direction.  Section 6 shows that even (7.6) implies generalized quadratic variance for symmetric log-concave measures.  Thus (7.5) identifies the correct modified defect, but it does not smuggle in a proof of the remaining thin-shell/KLS-strength statement.

## 8. Final obstruction map

The proposed proof route fails in two logically separate layers.

1. **Affine tangent failure.**  For every non-Gaussian convolution law, a linear-plus-polynomial perturbation has excess \(O(t)\) and squared Hoeffding defect \(O(t^2)\).  This disproves the stated inequalities outright, even for smooth strongly log-concave measures.

2. **Quadratic inverse problem.**  After projecting away affine functions or adding the exact Stein correction, pure quadratic modes demand
   \[
   \operatorname{Var}(X^TAX)\lesssim\operatorname{tr}(A^2),
   \]
   including the thin-shell conjecture at \(A=I\).

The four-copy residual measures nonlinear failure of additivity, but it is blind to the first variation of the Poincare excess along the affine eigenspace.  The affine--Stein term in (7.5) is the missing first-order defect.  Any viable convolution Lyapunov functional must track it, and any dimension-free inverse theorem on the remaining affine-orthogonal sector must still solve the generalized quadratic variance obstruction rather than assume it.
