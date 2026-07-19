# Linear \(H^{-1}\), minimal Stein fields, and the spectral-stability bottleneck

## 1. Verdict

Let \(\mu\) be an isotropic log-concave probability measure on
\(\mathbb R^n\), let \(a\in S^{n-1}\), and put

\[
 h_a(x)=a\cdot x,
 \qquad
 H_\mu(a)=
 \sup_g\frac{\left(\int h_a g\,d\mu\right)^2}
 {\int |\nabla g|^2\,d\mu}.
 \tag{1.1}
\]

Thus \(H_\mu(a)=\|h_a\|_{H^{-1}(\mu)}^2\); throughout, \(H_\mu(a)\)
denotes the *squared* negative-Sobolev norm.

The proposed dimension-free estimate

\[
 \boxed{H_\mu(a)\le C\quad\text{for every isotropic log-concave }
 \mu\text{ and every }a\in S^{n-1}}
 \tag{\mathrm{D}H^{-1}}
\]

is **not supplied by the Klartag--Fathi moment-map estimates**. It is an
operator-norm strengthening of the trace estimate used in the recent proof
of the thin-shell conjecture. I did not find a proof or a counterexample in
the primary sources checked as of July 16, 2026. The accurate status is
therefore that
\((\mathrm{D}H^{-1})\) is an unresolved operator strengthening of a now-known
average theorem.

There is a useful hierarchy:

\[
 \boxed{\mathrm{KLS}}
 \Longrightarrow
 \boxed{(\mathrm{D}H^{-1})}
 \Longrightarrow
 \boxed{\sum_{i=1}^n H_\mu(e_i)\le Cn}
 \Longrightarrow
 \boxed{\operatorname{Var}_\mu |X|^2\le Cn}.
 \tag{1.2}
\]

The first implication is immediate from \(H_\mu(a)\le C_P(\mu)\). The middle
trace estimate is Theorem 1.2 in Klartag--Lehec's
[Thin-shell bounds via parallel coupling](https://arxiv.org/abs/2507.15495).
The last implication is the log-concave \(H^{-1}\) inequality applied to
\(|x|^2-n\). None of the reverse implications in (1.2) is presently
available.

In particular, constructing a log-concave family with
\(H_{\mu_n}(a_n)\to\infty\) would give \(C_P(\mu_n)\to\infty\), hence would
disprove KLS. It is therefore unsurprising that no such family is known.

On the positive side, \((\mathrm{D}H^{-1})\) would close Target B of
[spectral_extremizer_variation.md](./spectral_extremizer_variation.md) in a
form *stronger* than requested: the KKT equations would not be needed. If
\(f\) is a normalized first eigenfunction and

\[
 \delta=\operatorname{dist}_{L^2(\mu)}(f,U)<1,
 \qquad U=\{a\cdot x:a\in\mathbb R^n\},
\]

then

\[
 \lambda_1(\mu)\ge \frac{1-\delta^2}{C}.
 \tag{1.3}
\]

Section 7 gives the exact two-line proof.

## 2. Exact Hilbert-space and Stein-field equivalence

Let

\[
 \mathcal G_\mu=
 \overline{\{\nabla g:g\in C_c^\infty(\mathbb R^n)\}}
 ^{L^2(\mu;\mathbb R^n)}.
\]

For a measure on a convex body, use smooth functions on the body and the
corresponding weak Neumann/no-flux convention. Constants are quotiented out.

### Proposition 2.1 (minimal Stein field)

For fixed \(a\), the following statements are equivalent:

1. \(H_\mu(a)<\infty\).
2. There is \(F_a\in L^2(\mu;\mathbb R^n)\) such that
   \[
    \int F_a\cdot\nabla g\,d\mu
       =\int (a\cdot x)g\,d\mu
       \quad\text{for every admissible }g.
    \tag{2.1}
   \]
3. In distributions, with the natural no-flux boundary condition when
   appropriate,
   \[
    -\operatorname{div}(\mu F_a)=(a\cdot x)\mu.
    \tag{2.2}
   \]

Moreover, there is a unique solution \(F_a^*\in\mathcal G_\mu\), and

\[
 \boxed{
 H_\mu(a)=\int |F_a^*|^2\,d\mu
 =\min_{F_a\text{ satisfying }(2.1)}\int|F_a|^2\,d\mu.}
 \tag{2.3}
\]

Consequently, the existence of *some* Stein field with
\(\int|F_a|^2d\mu\le C\) implies \(H_\mu(a)\le C\), and conversely
\(H_\mu(a)\le C\) produces the minimal Stein field with exactly that energy.

#### Proof

On the homogeneous Sobolev space define

\[
 \ell_a(g)=\int(a\cdot x)g\,d\mu.
\]

Its squared dual norm is exactly (1.1). If this is finite, Riesz
representation gives a unique \(F_a^*\in\mathcal G_\mu\) such that
\(\ell_a(g)=\langle F_a^*,\nabla g\rangle_{L^2(\mu)}\), and its squared norm
is \(H_\mu(a)\). Conversely, any field satisfying (2.1) bounds the functional
by Cauchy--Schwarz. Finally, if \(F_a\) is any solution, then
\(F_a-F_a^*\perp\mathcal G_\mu\), so

\[
 \|F_a\|_2^2=\|F_a^*\|_2^2+\|F_a-F_a^*\|_2^2.
\]

This proves (2.3), while (2.1) and (2.2) are the same weak identity.
\(\square\)

For a nondegenerate log-concave measure with finite covariance,
\(C_P(\mu)<\infty\). Hence

\[
 H_\mu(a)\le C_P(\mu)\int(a\cdot x)^2d\mu=C_P(\mu),
 \tag{2.4}
\]

so qualitative existence in Proposition 2.1 is not an issue; only the
dimension-free energy is.

### 2.1 Poisson representative

In the smooth setting, write

\[
 A=-\Delta+\nabla V\cdot\nabla
   =-\mu^{-1}\operatorname{div}(\mu\nabla)
\]

for the positive Langevin/Neumann generator. The weak solution of

\[
 Au_a=a\cdot x
 \tag{2.5}
\]

satisfies \(F_a^*=\nabla u_a\). Thus

\[
 H_\mu(a)=\langle h_a,A^{-1}h_a\rangle_{L^2(\mu)}.
 \tag{2.6}
\]

If \((\phi_k,\lambda_k)_{k\ge1}\) is a spectral resolution on the centered
subspace, then

\[
 H_\mu(a)=\sum_{k\ge1}
 \frac{|\langle h_a,\phi_k\rangle|^2}{\lambda_k}.
 \tag{2.7}
\]

This makes the nature of a possible obstruction transparent: a large value
requires a low eigenmode having appreciable correlation with a linear
functional.

### 2.2 Normalizations forced by isotropy

Testing (2.1) with \(g(x)=b\cdot x\) gives

\[
 \int F_a^*\,d\mu=a.
 \tag{2.8}
\]

Testing (1.1) with \(g=h_a\) gives \(H_\mu(a)\ge1\). Hence

\[
 H_\mu(a)
 =1+\int|F_a^*-a|^2d\mu.
 \tag{2.9}
\]

For smooth \(V\), the integrated Bochner identity applied to (2.5) yields

\[
 1=\int(Au_a)^2d\mu
 =\int\|D^2u_a\|_{\mathrm{HS}}^2d\mu
  +\int D^2V(F_a^*,F_a^*)d\mu,
 \tag{2.10}
\]

with an additional nonnegative Reilly boundary term on a convex Neumann
domain. Thus log-concavity controls the derivative of \(F_a^*\), but the
desired estimate asks to control the \(L^2\) oscillation in (2.9). Applying a
Poincare inequality here merely returns a dependence on \(C_P(\mu)\). This is
the precise point at which the elementary Poisson/Bochner proof stalls.

## 3. A common minimal Stein kernel and the operator formulation

The Riesz representative depends linearly on \(a\). Let \(u_i=u_{e_i}\), and
define the matrix field \(\tau_*\) by declaring its \(i\)-th row to be
\(\nabla u_i\). Then

\[
 F_a^*=\tau_*^T a.
 \tag{3.1}
\]

For every smooth vector field \(f=(f_1,\ldots,f_n)\),

\[
 \int x\cdot f\,d\mu
 =\int\langle\tau_*,\nabla f\rangle_{\mathrm{HS}}d\mu.
 \tag{3.2}
\]

Thus \(\tau_*\) is a Stein kernel, although it need not be symmetric or
positive semidefinite. Also

\[
 \mathbb E\tau_*=I_n.
 \tag{3.3}
\]

Define the deterministic positive-semidefinite matrix

\[
 K_\mu=\int\tau_*\tau_*^T d\mu.
 \tag{3.4}
\]

Then

\[
 \boxed{H_\mu(a)=a^TK_\mu a.}
 \tag{3.5}
\]

Consequently,

\[
 (\mathrm{D}H^{-1})
 \quad\Longleftrightarrow\quad
 K_\mu\preceq CI_n.
 \tag{3.6}
\]

This is the clean operator form of the question.

The kernel \(\tau_*\) is simultaneously minimal. If \(\tau\) is any other
\(L^2\) Stein kernel, then, for every \(a\),
\(\tau^Ta-\tau_*^Ta\) is orthogonal to all gradients. The two cross terms in
the expansion below therefore have zero quadratic form in every \(a\), and
hence vanish as a symmetric matrix. Thus

\[
 \int\tau\tau^T d\mu
 =K_\mu+\int(\tau-\tau_*)(\tau-\tau_*)^T d\mu
 \succeq K_\mu.
 \tag{3.7}
\]

In particular, \((\mathrm{D}H^{-1})\) is equivalent to the existence of a
single Stein kernel \(\tau\) satisfying

\[
 \int\tau\tau^Td\mu\preceq CI_n.
 \tag{3.8}
\]

The Hilbert-space construction then produces the optimal such kernel.

## 4. What is now known: trace control, not operator control

Klartag--Lehec prove that for every isotropic log-concave \(\mu\),

\[
 \sum_{i=1}^n\|x_i\|_{H^{-1}(\mu)}^2\le Cn.
 \tag{4.1}
\]

In the present notation this is exactly

\[
 \boxed{\operatorname{tr}K_\mu
       =\sum_{i=1}^nH_\mu(e_i)
       =\int\|\tau_*\|_{\mathrm{HS}}^2d\mu
       \le Cn.}
 \tag{4.2}
\]

This is optimal in order, since (3.3) and Jensen give
\(\int\|\tau_*\|_{\mathrm{HS}}^2d\mu\ge n\). It also says that for a uniformly
random direction \(A\in S^{n-1}\),

\[
 \mathbb E_A H_\mu(A)=\frac{\operatorname{tr}K_\mu}{n}\le C.
 \tag{4.3}
\]

It does **not** bound a prescribed or adversarial direction: a matrix with
trace \(Cn\) may have one eigenvalue of order \(n\). The direction supplied by
the tensor-rigidity argument in
[spectral_extremizer_variation.md](./spectral_extremizer_variation.md) is
selected by the eigenfunction, so (4.3) cannot be substituted for (3.6).

The usual \(H^{-1}\) inequality (Proposition 10 in Barthe--Klartag,
[Spectral gaps, symmetries and log-concave perturbations](https://arxiv.org/abs/1907.01823))
gives, for centered \(f\) whose partial derivatives are centered,

\[
 \operatorname{Var}_\mu f
 \le\sum_i\|\partial_i f\|_{H^{-1}(\mu)}^2.
 \tag{4.4}
\]

Taking \(f(x)=|x|^2-n\) yields

\[
 \operatorname{Var}_\mu |X|^2
 \le4\sum_iH_\mu(e_i).
 \tag{4.5}
\]

Notice the direction: a thin-shell estimate by itself gives no upper bound on
the right-hand side. Klartag--Lehec prove the stronger trace bound (4.1) on
their way to thin shell; thin shell cannot simply be run backwards to prove
directional control.

## 5. Moment-map kernel: exact benefit and exact limitation

Under the regularity assumptions in Fathi's
[Stein kernels and moment maps](https://arxiv.org/abs/1804.04699), let
\(\varphi\) be the moment potential, so that \(\nabla\varphi\) pushes
\(e^{-\varphi(y)}dy\) to \(\mu\). Then

\[
 \tau_{\mathrm{mm}}(x)
 =D^2\varphi(\nabla\varphi^*(x))
 \tag{5.1}
\]

is a symmetric positive-semidefinite Stein kernel. Therefore

\[
 F_a(x)=\tau_{\mathrm{mm}}(x)a
 \tag{5.2}
\]

is an admissible field in Proposition 2.1, and

\[
 H_\mu(a)
 \le\int|\tau_{\mathrm{mm}}a|^2d\mu
 =\int a^T\tau_{\mathrm{mm}}^2a\,d\mu.
 \tag{5.3}
\]

For isotropic log-concave \(\mu\), the Klartag estimate quoted as Fathi's
Proposition 3.2 gives, for \(p\ge1\),

\[
 \left\|a^T\tau_{\mathrm{mm}}a\right\|_{L^p(\mu)}
 \le 8p^2.
 \tag{5.4}
\]

But the required quantity in (5.3) is

\[
 |\tau_{\mathrm{mm}}a|^2
 =a^T\tau_{\mathrm{mm}}^2a,
 \tag{5.5}
\]

not \((a^T\tau_{\mathrm{mm}}a)^2\). Even for a symmetric PSD matrix, these
are different quantities, with

\[
 (a^T\tau a)^2\le a^T\tau^2a.
 \tag{5.6}
\]

Thus (5.4) bounds the smaller expression.

### 5.1 Algebraic audit: diagonal moments cannot control a column

Let \(U\) be uniform on \(S^{n-1}\) and set

\[
 T=nUU^T.
 \tag{5.7}
\]

Then \(T\succeq0\), \(T=T^T\), and \(\mathbb ET=I_n\). For every fixed unit
\(a\),

\[
 a^TTa=n(a\cdot U)^2,
 \qquad
 \|a^TTa\|_{L^p}\le 2p,
 \tag{5.8}
\]

while

\[
 \mathbb E|Ta|^2
 =n^2\mathbb E(a\cdot U)^2=n.
 \tag{5.9}
\]

This is not asserted to be the Stein kernel of a log-concave measure. Its role
is narrower and decisive: symmetry, positivity, mean identity, and uniform
directional quadratic-form moments do not logically imply a column-norm
bound. Some additional moment-map PDE structure would have to be used.

### 5.2 The weighted Poincare inequality does not repair the gap

The moment-map kernel satisfies

\[
 \operatorname{Var}_\mu g
 \le\int\nabla g^T\tau_{\mathrm{mm}}\nabla g\,d\mu.
 \tag{5.10}
\]

The Stein identity and pointwise PSD Cauchy--Schwarz give

\[
 |\mathbb E[h_ag]|^2
 \le
 \mathbb E[a^T\tau_{\mathrm{mm}}a]\,
 \mathbb E[\nabla g^T\tau_{\mathrm{mm}}\nabla g]
 =\mathbb E[\nabla g^T\tau_{\mathrm{mm}}\nabla g].
 \tag{5.11}
\]

The last energy is weighted and can be larger than
\(\mathbb E|\nabla g|^2\). Inequality (5.10) has the same direction and gives
no comparison with the unweighted Dirichlet energy. Hence this route also
stops short of (1.1).

Finally, failure to bound the moment-map column would not itself disprove
\((\mathrm{D}H^{-1})\): by (3.7), the minimal field is the gradient projection
of the moment-map field and may have much smaller energy.

## 6. Sanity checks on canonical log-concave families

No growth occurs in the standard extremal test families.

### 6.1 Gaussian and products

For the standard Gaussian, \(u_a(x)=a\cdot x\), hence \(F_a^*=a\) and

\[
 H_{\gamma_n}(a)=1.
\]

For a product of centered one-dimensional laws
\(\mu=\mu_1\otimes\cdots\otimes\mu_n\), let \(A_i u_i=x_i\). Then
\(u_a=\sum_i a_i u_i(x_i)\), and

\[
 H_\mu(a)=\sum_i a_i^2H_{\mu_i}(1).
 \tag{6.1}
\]

One-dimensional isotropic log-concavity has a universal Poincare constant,
so (6.1) is dimension-free. For the isotropic cube
\([-\sqrt3,\sqrt3]^n\), the one-dimensional Stein kernel is
\(\tau(x)=(3-x^2)/2\), giving the exact value

\[
 H_\mu(a)=\mathbb E\tau(X)^2=\frac65.
 \tag{6.2}
\]

### 6.2 The isotropic simplex

Let \(N=n+1\), let \(U\) be uniform on the probability simplex
\(\{u_i\ge0:\sum_i u_i=1\}\), and work in
\(E=\{z\in\mathbb R^N:\sum_i z_i=0\}\). Then

\[
 X=\sqrt{N(N+1)}\left(U-\frac1N\mathbf1\right)
\]

is isotropic in the \(n\)-dimensional space \(E\). The Dirichlet integration
by parts formula gives the symmetric PSD Stein kernel

\[
 \tau(X)=(N+1)(\operatorname{diag}U-UU^T)|_E.
 \tag{6.3}
\]

For every \(a\in E\) with \(|a|=1\), a direct fourth-moment calculation gives

\[
 \mathbb E|\tau a|^2=\frac{2(N+1)}{N+3}<2.
 \tag{6.4}
\]

Thus \(H_\mu(a)<2\) for every direction of the simplex. This also illustrates
the relevant distinction: here the full column happens to be controlled, but
that conclusion uses the explicit kernel, not only moments of \(a^T\tau a\).

More generally, if the orthogonal symmetry group of \(\mu\) acts irreducibly,
then uniqueness of the minimal kernel forces \(K_\mu\) to be scalar. The trace
bound (4.2) then implies \((\mathrm{D}H^{-1})\) for that measure. A potential
counterexample must therefore be substantially anisotropic beyond covariance.

## 7. Exact closure of Target B

Let \(f\) be a centered, \(L^2(\mu)\)-normalized first eigenfunction:

\[
 \int f^2d\mu=1,
 \qquad
 \int|\nabla f|^2d\mu=\lambda_1(\mu)=:\lambda.
 \tag{7.1}
\]

Because \(\mu\) is isotropic, the map \(a\mapsto a\cdot x\) is an isometry
from \(\mathbb R^n\) onto the linear subspace \(U\subset L^2(\mu)\). Write

\[
 P_Uf=a\cdot x,
 \qquad
 \delta=\|f-P_Uf\|_2.
\]

Orthogonality gives

\[
 |a|^2=\|P_Uf\|_2^2=1-\delta^2.
 \tag{7.2}
\]

If \(\delta<1\), set \(b=a/|a|\). Then

\[
 \left(\int f(x)(b\cdot x)d\mu(x)\right)^2
 =|a|^2=1-\delta^2.
 \tag{7.3}
\]

Use \(g=f\) in the definition of \(H_\mu(b)\):

\[
 1-\delta^2
 \le H_\mu(b)\int|\nabla f|^2d\mu
 =H_\mu(b)\lambda.
 \tag{7.4}
\]

Under \((\mathrm{D}H^{-1})\), this proves

\[
 \boxed{\lambda\ge\frac{1-\delta^2}{C}.}
 \tag{7.5}
\]

Equivalently, using a Stein field of energy at most \(C\), (7.3) and (2.1)
give

\[
 \sqrt{1-\delta^2}
 =\left|\int F_b\cdot\nabla f\,d\mu\right|
 \le\sqrt C\sqrt\lambda.
\]

Therefore Target B can take, for example,

\[
 \varepsilon_0=\frac12,
 \qquad
 c_0=\frac{3}{4C}.
 \tag{7.6}
\]

No normalized KKT equation is used. Mere \(L^2\)-closeness plus the
Rayleigh energy of the eigenfunction suffices once the directional
\(H^{-1}\) estimate is available.

With only the trace theorem (4.2), one gets at best
\(H_\mu(b)\le Cn\) for an adversarial \(b\), and hence only
\(\lambda\gtrsim (1-\delta^2)/n\). That does not close Target B.

## 8. Bottom line for the spectral-extremizer program

The linear \(H^{-1}\) proposal is mathematically exact and would completely
remove Target B as a bottleneck. It is equivalent to any of the following
uniform operator statements:

\[
 \sup_{|a|=1}H_\mu(a)\le C,
 \qquad
 \int|F_a^*|^2d\mu\le C,
 \qquad
 K_\mu\preceq CI,
\]

or to the existence of a common Stein kernel with
\(\int\tau\tau^Td\mu\preceq CI\).

What is currently proved is the Hilbert--Schmidt/trace analogue

\[
 \int\|\tau_*\|_{\mathrm{HS}}^2d\mu\le Cn.
\]

The moment-map estimates control \(a^T\tau a\), while Target B needs a bound
on \(|\tau a|\), or more optimally on the gradient projection of that column.
Conflating these quantities is precisely the invalid step. Proving the
operator upgrade would be a genuine new result intermediate between KLS and
thin shell; producing a growing example would disprove KLS.
