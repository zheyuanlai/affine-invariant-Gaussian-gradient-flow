# Gradient multiplicity and rigidity in the Gaussian observation model

## 0. Verdict

There is a useful exact near-equality decomposition behind the convolution
Poincare estimate, but **multiplicity does not by itself force a Gaussian or
product factor**.  The decomposition separates three different defects.  One
of them is a Gaussian-noise defect and has a dimension-free Hermite stability
theorem; another is a Poincare defect for the original measure and has no
known factor-rigidity theorem.  Identifying its high-multiplicity near
eigenspace with a product factor would be a new spectral-rigidity hypothesis,
not a consequence of tensorization.

There are two concrete warnings.

1. Irreducible radial log-concave measures have large spectral
   multiplicities for representation-theoretic reasons, without possessing
   a nontrivial product factor.
2. For a median Euclidean ball under a Gaussian measure, the Fisher matrix is
   a scalar matrix and the wedge-Poincare bound is asymptotically paid by the
   angular Hessian with an explicitly small excess.  An arbitrarily small
   radial quartic perturbation removes every exact Gaussian/product factor
   while preserving all of these inequalities to arbitrary prescribed
   accuracy in each fixed dimension.

The second example does not rule out a *quantitative approximate-factor
theorem at fixed positive* `r_0` and fixed `alpha=s/K`: its excess and Fisher
trace obey an explicit tradeoff.  It does rule out the proposed qualitative
conclusion of a genuine factor and shows exactly which parameter dependence a
valid inverse theorem would have to retain.

The companion result `general_angular_stability.md` was independently
audited.  Its source-halfspace approximation, dimension-free
Hilbert--Schmidt perturbation, projection change, and `t^{-2}` rescaling are
all valid.  Section 5 below records the audit and the exact way in which that
result controls the angular, but not the longitudinal, Hessian.

---

## 1. Exact defect decomposition for convolution Poincare

Let `mu` be isotropic and let `C_P(mu)=K<infinity`.  Let

\[
 X\sim\mu,\qquad G\sim\gamma_n,\qquad Y=X+\sqrt{s}\,G,
 \qquad q_s={\cal L}(Y),                                      \tag{1.1}
\]

where `X` and `G` are independent.  For a scalar
`phi in W^{1,2}(q_s)`, define

\[
 \Phi(x)=\mathbb E_G\phi(x+\sqrt{s}G),\qquad
 J(x)=\mathbb E_G\nabla\phi(x+\sqrt{s}G)=\nabla\Phi(x).       \tag{1.2}
\]

Put

\[
\begin{aligned}
 {cal D}_G(\phi)
 &=s\mathbb E|\nabla\phi(Y)|^2
   -\mathbb E_X\operatorname {Var}_G
        \big(\phi(X+\sqrt{s}G)\big),\\
 {cal D}_\mu(\phi)
 &=K\mathbb E_\mu|J|^2-\operatorname {Var}_\mu(\Phi),\\
 {cal D}_J(\phi)
 &=\mathbb E|\nabla\phi(Y)|^2-\mathbb E_\mu|J|^2.
                                                               \tag{1.3}
\end{aligned}
\]

Gaussian Poincare, Poincare for `mu`, and conditional Jensen show that all
three terms are nonnegative.  The variance decomposition gives the exact
identity

\[
\boxed{
 (K+s)\mathbb E|\nabla\phi|^2-\operatorname {Var}_{q_s}(\phi)
 ={cal D}_G(\phi)+{cal D}_\mu(\phi)+K{cal D}_J(\phi).}
                                                               \tag{1.4}
\]

Thus near equality in the tensorized bound
`C_P(q_s)<=K+s` forces three simultaneous statements, not one:

* near equality in Gaussian Poincare on most conditional translates;
* near equality in the original `mu`-Poincare inequality for `Phi`;
* near equality in conditional Jensen for the gradient.

The first item has an exact dimension-free inverse.  If
`f(G)=sum_{d>=0}f_d(G)` is its Hermite decomposition, then

\[
 s\mathbb E|\nabla f|^2-\operatorname {Var}(f)
 =\sum_{d\ge2}(d-1)\|f_d\|_2^2.                       \tag{1.5}
\]

Consequently a small `D_G` makes the conditional functions close in `L^2`
to affine functions of `G`.  It does not make the affine slopes attached to
distant values of `X` agree.  The term `D_J` aligns a slope only across the
Gaussian translate based at a fixed `X`; when different translates have
negligible overlap, a separate gluing argument is still needed.  Finally,
small `D_mu` asks for near extremizers of the original Poincare inequality.
Classifying a high-dimensional space of such extremizers as a product factor
is not a standard equality theorem for arbitrary log-concave measures.

### 1.1 Exact equality

The equality case is elementary and is worth separating from near equality.
If `D_J(phi)=0`, then, for `mu`-almost every `x`,

\[
 \nabla\phi(x+\sqrt{s}G)=J(x)\quad\text{for Gaussian-a.e. }G.  \tag{1.6}
\]

Every Gaussian translate has a strictly positive Lebesgue density.  Comparing
two admissible values of `x` therefore shows that `J(x)` is a single constant
and that `phi` is affine Lebesgue-a.e.  For
`phi(y)=theta\mathbin\cdot y+c`, isotropy gives

\[
 {cal D}_\mu(\phi)=(K-1)|\theta|^2.                           \tag{1.7}
\]

Hence a nonconstant exact extremizer of the particular upper bound `K+s`
can occur only when `K=1`, and it is affine.  This exact statement does not
have a robust product conclusion: robustness requires quantitative gluing of
the conditional affine approximations and quantitative classification of
the `mu`-near eigenspace.

---

## 2. Applying the defect decomposition to all wedge components

Let `F:R^n->R^n` and let `Y'` be an independent copy of `Y`.  For fixed
`a=F(Y')`, apply (1.4) componentwise to the exterior-product-valued map

\[
             \phi_a(y)=a\wedge F(y).                           \tag{2.1}
\]

Summing the scalar decompositions over an orthonormal basis of
`Lambda^2 R^n` and then averaging over `Y'` gives an exact decomposition of
the deficit in the wedge-Poincare step.  In particular, if that step is
nearly sharp, each of the three nonnegative defects in (1.4), averaged over
all random wedge components, is small.

There is a clean exact obstruction to simultaneous equality.  Suppose the
linear span of the essential range of `F` has dimension at least two and
every map `a wedge F` appearing in (2.1) is affine.  Taking two weak
derivatives gives

\[
             a\wedge \partial_{ij}F=0                           \tag{2.2}
\]

for almost every `a` in the range.  Two nonparallel such vectors force
`partial_ij F=0`; hence `F` is affine.  In the posterior application

\[
       F=\sqrt{s}\,\nabla h,\qquad 0\le h\le\pi.                \tag{2.3}
\]

Since `q_s` has a positive density on all of `R^n`, an affine `F` makes `h`
a globally bounded quadratic polynomial.  It must therefore be constant,
so `F=0`.  Thus exact nonzero high-rank simultaneous equality is impossible;
it does not produce a nontrivial factor.

This algebra is not quantitatively stable without more information.  An
`L^2(q_s)`-small error permits different almost-affine descriptions on
regions connected only through low-amplitude zones.  That is precisely the
phase-cell obstruction.  Replacing this missing quantitative statement by
"many near eigenfunctions imply a Gaussian factor" would assume the desired
spectral rigidity rather than prove it.

Large multiplicity alone is especially weak evidence for a factor.  A
strictly convex radial potential

\[
 V_\epsilon(x)=\frac{|x|^2}{2}+\frac{\epsilon}{4n}|x|^4,
 \qquad \epsilon>0,                                           \tag{2.4}
\]

defines a full-dimensional, smooth, strongly log-concave measure which has
no nontrivial orthogonal product decomposition.  Nevertheless its diffusion
commutes with `O(n)`, so every eigenfunction in the first spherical-harmonic
sector occurs with multiplicity `n`.  For sufficiently small `epsilon` this
sector is the perturbation of the `n`-dimensional first Gaussian eigenspace.
Thus even exact spectral multiplicity of the original diffusion does not
encode product structure.

---

## 3. The exact `alpha` bookkeeping

Use the posterior notation from `gradient_wedge_poincare.md`:

\[
 g(y)=\mathbb P(B=1\mid Y=y),\quad
 h=2\arcsin\sqrt g,\quad F=\sqrt{s}\,\nabla h,\quad
 R=\mathbb E[FF^{\mathsf T}].                                 \tag{3.1}
\]

Assume

\[
 \operatorname {tr}R\ge r_0,\qquad R\preceq s^{-1}I,
 \qquad s=\alpha K.                                           \tag{3.2}
\]

Then `kappa<=1/(sr_0)`.  Once `sr_0>2`, the wedge inequality gives

\[
\boxed{
 \mathbb E\|\nabla^2h\|_{\mathrm{HS}}^2
 \ge
 \left(1-\frac{2}{sr_0}\right)
 \frac{r_0}{s(K+s)}
 =\left(1-\frac{2}{\alpha Kr_0}\right)
 \frac{r_0}{\alpha(1+\alpha)K^2}.}                            \tag{3.3}
\]

The harmless factor in parentheses is retained here because dropping it
can conceal dependence on `alpha` and `r_0`.

Write `z=Phi^{-1}(g)` and

\[
 a(z)=\frac{\varphi(z)}{\sqrt{\Phi(z)(1-\Phi(z))}}.
\]

At points where `nabla z` is nonzero, put
`e=nabla z/|nabla z|` and `P=I-ee^T`.  If `H=nabla^2z` and
`w=|nabla z|^2`, then the three orthogonal blocks of `nabla^2h` give the
exact identity

\[
\boxed{
\begin{aligned}
 \|\nabla^2h\|_{\mathrm{HS}}^2
 &=a^2\|PHP\|_{\mathrm{HS}}^2
   +2a^2|PHe|^2
   +\big(a\,e^{\mathsf T}He+a'w\big)^2.
                                                               \tag{3.4}
\end{aligned}}
\]

The first two terms are angular (transverse-transverse and mixed).  The last
term is longitudinal.  It is important not to call the whole
`a nabla^2 z` term angular: its `e-e` entry can cancel or reinforce the
amplitude derivative `a'w`.

The posterior is `s^{-1}`-strongly log-concave, and the sharp centroid
inequality gives `w<=s^{-1}`.  Since `a'` is bounded,

\[
 \mathbb E\|a'(z)\nabla z\nabla z^{\mathsf T}\|_{\mathrm{HS}}^2
 \le \frac{C}{s^2}=\frac{C}{\alpha^2K^2}.                     \tag{3.5}
\]

Relative to the lower scale in (3.3), the available pure-amplitude budget is

\[
 \frac{C/(\alpha^2K^2)}{r_0/[\alpha(1+\alpha)K^2]}
 =\frac{C(1+\alpha)}{\alpha r_0}.                              \tag{3.6}
\]

Thus at small `alpha` the longitudinal budget is *larger*, not smaller,
than the wedge lower bound by a factor of order `(alpha r_0)^{-1}`.  No
angular conclusion follows from the Bernstein-size estimate alone.

---

## 4. A radial high-rank near-saturation model

This section gives a model in which all Fisher directions coexist because
of rotational symmetry rather than because of independent factors.

Let `mu=gamma_n`, so `K=1`, fix `s>0`, and put `beta=1+s`.  Let

\[
 B_n=1_{\{|X|^2\le m_n\}},                                    \tag{4.1}
\]

where `m_n` is a median of `chi_n^2`.  Then `g=1/2`, and the posterior
function `g_n(y)` and `h_n(y)=2arcsin sqrt(g_n(y))` are radial.  Hence

\[
 R_n=\mathbb E[(\sqrt{s}\nabla h_n)
                (\sqrt{s}\nabla h_n)^{\mathsf T}]
     =\frac{r_n(s)}{n}I,\qquad \mathbb E F=0.                 \tag{4.2}
\]

Define

\[
 U_n=\frac{|X|^2-n}{\sqrt{2n}},\qquad
 V_n=\frac{|Y|^2-\beta n}{\beta\sqrt{2n}}.                    \tag{4.3}
\]

The joint central limit theorem gives

\[
 (U_n,V_n)\Longrightarrow (U,V),\qquad
 \operatorname {Corr}(U,V)=\frac1\beta.                       \tag{4.4}
\]

Since `(m_n-n)/sqrt(2n)->0`, the conditional posterior has the local-CLT
limit

\[
 g_n(Y)\Longrightarrow \Phi(-\lambda V),\qquad
 \lambda=\frac{1}{\sqrt{\beta^2-1}}
        =\frac1{\sqrt{s(s+2)}}.                                \tag{4.5}
\]

The same local CLT differentiated twice (or the explicit noncentral
`chi^2` density) gives convergence of the radial first and second derivative
energies below.  With

\[
 a(z)=\frac{\varphi(z)}{\sqrt{\Phi(z)(1-\Phi(z))}},\qquad
 Z\sim N(0,1),                                                 \tag{4.6}
\]

one obtains

\[
\boxed{
 r_n(s)\longrightarrow
 r_\infty(s)=\frac{2s\lambda^2}{\beta}
                 \mathbb E a(\lambda Z)^2
             =\frac{2}{\beta(s+2)}\mathbb E a(\lambda Z)^2.} \tag{4.7}
\]

In particular `r_infty(s)>0` for every fixed `s`, and (4.2) obeys
`R_n preceq s^{-1}I` for all sufficiently large `n`.

For a radial function, the Hessian eigenvalues are `h_rr` once and
`h_r/r` with multiplicity `n-1`.  The two limiting contributions are

\[
\begin{aligned}
 A_s&:=\lim_{n\to\infty}
  \mathbb E (n-1)\left(\frac{\partial_rh_n}{r}\right)^2
   =\frac{2\lambda^2}{\beta^2}\mathbb E a(\lambda Z)^2
   =\frac{r_\infty(s)}{s(1+s)},\\
 L_s&:=\lim_{n\to\infty}
  \mathbb E(\partial_{rr}h_n)^2
   =\frac{4\lambda^4}{\beta^2}
      \mathbb E a'(\lambda Z)^2.                              \tag{4.8}
\end{aligned}
\]

Here `A_s` is purely angular and `L_s` is longitudinal.  On the other hand,
because `q_s=N(0,beta I)`, `R_n=(r_n/n)I`, and the mean of `F` vanishes,
the exact wedge-Poincare lower bound is

\[
 \mathbb E\|\nabla^2h_n\|_{\mathrm{HS}}^2
 \ge\frac{r_n(s)}{s(1+s)}.                                   \tag{4.9}
\]

Thus the angular term alone converges exactly to the lower bound.  The
relative excess is

\[
 \frac{L_s}{A_s}
 =2\lambda^2\frac{\mathbb E a'(\lambda Z)^2}
                       {\mathbb E a(\lambda Z)^2}.             \tag{4.10}
\]

The function `a` is even and smooth at zero, so `a'(0)=0`.  Dominated
convergence with the Gaussian tails yields

\[
 \frac{L_s}{A_s}=O(s^{-4}),\qquad
 r_\infty(s)\sim\frac{4}{\pi s^2}\quad(s\to\infty).        \tag{4.11}
\]

This is a simultaneous high-rank near-saturation mechanism with no phase
decomposition: all directions are the angular first harmonics of one radial
field.  It also records the unavoidable parameter tradeoff.  At
`alpha=s`, an excess `O(alpha^{-4})` is obtained with
`r_0 asymp alpha^{-2}`.  Therefore the example refutes parameter-free
qualitative factor rigidity, but it does not refute a stability gap whose
modulus is allowed to depend on a fixed lower bound for `r_0` and on a fixed
`alpha`.

### 4.1 Removing every exact factor

For each `n`, perturb the Gaussian to

\[
 d\mu_{n,\epsilon}(x)=Z_{n,\epsilon}^{-1}
 \exp\left(-\frac{|x|^2}{2}-\frac{\epsilon}{4n}|x|^4\right)dx, \tag{4.12}
\]

and apply the scalar radial dilation which makes the measure isotropic.
The resulting measure is smooth, full-dimensional, and strongly
log-concave.  It has no nontrivial product decomposition under any
orthogonal splitting: for `R^n=E oplus E^perp`, the quartic term contains

\[
 \frac{\epsilon}{2n}|x_E|^2|x_{E^\perp}|^2,                   \tag{4.13}
\]

which cannot be written as a sum of an `E`-function and an
`E^perp`-function.  In fact no invertible affine change produces a product:
after such a change the quartic term is a nonzero multiple of
`(z^TQz)^2` for a positive definite `Q`, whose mixed fourth derivatives
cannot vanish across a nontrivial splitting.  In particular the measure has
no Gaussian factor.

For fixed `n` and `s>0`, as `epsilon downarrow0`, the measures, their
isotropic dilation factors, their median radial balls, and the Gaussian
convolutions converge to their Gaussian counterparts.  Differentiating the
Gaussian convolution kernel up to order two and truncating first to a large
ball gives convergence of `r_n(s)` and of both Hessian energies.  The tails
are uniform because these perturbations are uniformly strongly log-concave.
Consequently, for any prescribed sequence of tolerances, one can choose a
positive diagonal sequence `epsilon_n downarrow0` so that the irreducible
measures (4.12) reproduce (4.7)--(4.10) to those tolerances.

This supplies full-dimensional irreducible countermodels to the assertion
that near simultaneous saturation forces a **genuine** factor.  What
survives is only approximate Gaussianity of this particular example.

---

## 5. Independent audit of `general_angular_stability.md`

The theorem in that file is valid, conditional only on the already proved
clean-room spatial stability theorem and the halfspace-pullback contraction
theorem.  The four possible failure points check as follows.

### 5.1 Source halfspace approximation

After standardizing to strong-convexity parameter one, set

\[
 Y=\langle u,T(G)\rangle,\qquad Z=\langle u,G\rangle,qquad
 \mathbb E|Y-Z|^2\le\zeta.                                    \tag{5.1}
\]

If `E={Y>=a}` and `F={Z>=b}` have the same mass, compare first with
`F_a={Z>=a}`.  For every `h>0`,

\[
 \mathbb P(E\mathbin\triangle F_a)
 \le\frac{\zeta}{h^2}+\frac{2h}{\sqrt{2\pi}}.                 \tag{5.2}
\]

The Gaussian threshold sets `F_a,F` are nested and equal mass of `E,F`
implies

\[
 \mathbb P(F_a\mathbin\triangle F)
 \le\mathbb P(E\mathbin\triangle F_a).                        \tag{5.3}
\]

Choosing `h=zeta^(1/3)` proves the claimed `C zeta^(1/3)` error.  Combining
it with the target-space symmetric-difference estimate gives exactly

\[
 \rho=C_\delta\left\{\sqrt\varepsilon+
  (\varepsilon\sqrt{\log(e/\varepsilon)})^{1/3}\right\}.       \tag{5.4}
\]

Only one-dimensional Gaussian anti-concentration is used.

### 5.2 Dimension-free Hilbert--Schmidt perturbation

Let `sigma=1_A-1_B`.  Equal masses give `E sigma=0` and
`E sigma^2=rho`.  For every symmetric `M` with `||M||HS=1`, one-strong
log-concavity and Poincare give

\[
 \operatorname {Var}(X^{\mathsf T}MX)
 \le4\mathbb E|MX|^2
 =4\operatorname {tr}(M^2\operatorname {Cov}X)\le4.            \tag{5.5}
\]

Therefore

\[
 |\operatorname {tr}M(D_A-D_B)|
 \le2\sqrt\rho.                                               \tag{5.6}
\]

Hilbert--Schmidt duality over symmetric matrices gives
`||D_A-D_B||HS<=2sqrt(rho)`.  This step is genuinely dimension free; it
pairs `Cov X preceq I` with `tr M^2=1`, not with `tr Cov X`.  The source
indicator `B` need not be measurable with respect to `X`: Poincare is
applied only to the quadratic function of `X`, followed by Cauchy--Schwarz
on the joint source space.

Taking `M=theta theta^T` and using
`E(1_B-g)^2=g(1-g)<=1/4` similarly gives `||D_B||op<=1`.

### 5.3 Projection change

Centrality and `|v_A-v_B|<=sqrt(rho)` imply
`|u_A-u_B|<=C_delta sqrt(rho)`.  Since
`P_A-P_B=u_Bu_B^T-u_Au_A^T` has rank at most two,

\[
 \|(P_A-P_B)D_B\|_{\mathrm{HS}}
 \le\|P_A-P_B\|_{\mathrm{HS}}\|D_B\|_{\mathrm{op}}
 \le C_\delta\sqrt\rho.                                      \tag{5.7}
\]

There is no hidden `sqrt(n)`.

### 5.4 Scaling

For a `t`-strongly log-concave target, standardize by
`X_tilde=sqrt(t)(X-m)`.  Then

\[
 \widetilde D=tD,\qquad \widetilde v=\sqrt t\,v,\qquad
 \widetilde P=P,\qquad \widetilde\varepsilon=\varepsilon.     \tag{5.8}
\]

Thus `||P Dtilde||HS^2<=C Omega` becomes

\[
 \boxed{\|PD\|_{\mathrm{HS}}^2
 \le C_\delta t^{-2}\Omega_\delta(\varepsilon).}             \tag{5.9}
\]

The statement remains intrinsic on a lower-dimensional affine hull by
performing the Gaussian transport and all matrix norms on its parallel
linear space.

### 5.5 Translation to the heat Hessian

At heat variance `s`, the posterior is `t=1/s` strongly log-concave.  If
`D` is its centered set-covariance and `I(g)=varphi(z)`, the exact heat
identity is

\[
 P\nabla_y^2z=\frac{PD}{s^2I(g)}.                              \tag{5.10}
\]

For central `g`, (5.9) therefore gives

\[
 \|P\nabla^2z\|_{\mathrm{HS}}^2
 \le\frac{C_\delta}{s^2}\Omega_\delta(\varepsilon),
 \qquad
 \|P\nabla^2h\|_{\mathrm{HS}}^2
 \le\frac{C_\delta}{s^2}\Omega_\delta(\varepsilon).         \tag{5.11}
\]

Here the left projection annihilates the amplitude term in (3.4).  At
`s=alpha K`, the ratio of this angular upper scale to (3.3) is

\[
 \frac{C\Omega/(\alpha^2K^2)}
      {r_0/[\alpha(1+\alpha)K^2]}
 =\frac{C(1+\alpha)}{\alpha r_0}\Omega_\delta(\varepsilon).   \tag{5.12}
\]

Thus angular energy is negligible compared with the wedge lower bound only
when

\[
             \Omega_\delta(\varepsilon)\ll
             \frac{\alpha r_0}{1+\alpha}.                      \tag{5.13}
\]

Even in that regime, (3.5)--(3.6) leave enough longitudinal amplitude
budget to pay the lower bound, especially for small `alpha`.  The angular
stability theorem is correct and useful, but it does not supply the missing
longitudinal exclusion or the gluing of phase cells.

---

## 6. What a viable rigidity lemma would have to prove

A noncircular inverse theorem can still be useful, but it must retain all
of the following quantitative data.

1. It must fix lower bounds on `r_0` and on the relevant central posterior
   mass, and retain explicit dependence on `alpha=s/K`.
2. It must turn the Gaussian conditional affine approximations from
   `D_G` into a single global approximation.  This is an overlap/gluing
   problem, not a consequence of Hermite stability on each translate.
3. It must control the longitudinal square
   `(a e^THe+a'w)^2`, including possible cancellation, rather than merely
   bound `a'^2w^2`.
4. It must rule out rotational high-multiplicity mechanisms such as Section
   4, not merely irreducible phase cells.
5. Any factor conclusion must be approximate and metrized.  Exact product
   structure is destroyed by arbitrarily small strictly convex couplings.
6. It may not assume that a high-multiplicity near eigenspace of an arbitrary
   log-concave diffusion splits the measure.  That assertion is neither a
   general Obata theorem nor a known consequence of log-concavity and is at
   least as strong as the missing part of this route.

The sharp current conclusion is therefore a trichotomy: substantial angular
energy, substantial longitudinal amplitude energy, or a small aggregate
convolution-Poincare defect requiring a new quantitative gluing theorem.
None of the three existing estimates forces a genuine Gaussian/product
factor.
