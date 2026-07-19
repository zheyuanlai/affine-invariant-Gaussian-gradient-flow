# Spectral rigidity audit: the natural one-form estimate is KLS-equivalent

This note audits the following tempting route.  If a normalized first
eigenfunction has a very small eigenvalue, the Bochner identity says that its
gradient is an almost-parallel curl-free vector field.  One might try to rule
out such a field by a dimension-free rigidity theorem.  The calculation below
shows that the natural rigidity theorem is quantitatively equivalent to the
original scalar Poincare problem.  It is therefore not an independent lemma.

## 1. Formulation that also covers non-attainment

Let `mu` be an isotropic log-concave probability on `R^n` with full-dimensional
support.  Let

```
q(f,g) = int <grad f,grad g> dmu
```

be its closed Dirichlet form and let `A=-L` be the associated nonnegative
self-adjoint operator on `L^2(mu)`.  Constants form the kernel.  Write

```
lambda = inf spectrum(A restricted to 1^perp),
C_P(mu) = 1/lambda.
```

The coordinate functions belong to the form domain, have variance one, and
have Dirichlet energy one.  Hence `0 < lambda <= 1` whenever `C_P(mu)` is
finite.  (For a log-concave probability the spectral gap is positive in each
fixed finite dimension; no dimension-free lower bound is being used here.)

For `f in Dom(A)` define

```
m_f = int grad f dmu,
N(f) = int |grad f-m_f|^2 dmu.
```

Let `R_mu` be the optimal constant in

```
N(f) <= R_mu ||Af||_2^2                         (1)
```

for all mean-zero `f in Dom(A)`.  Then

```
C_P(mu)-1 <= R_mu <= C_P(mu).                  (2)
```

All constants in (2) are exact and independent of the dimension.

### Proof of the upper bound

The spectral theorem gives, on `1^perp`,

```
int |grad f|^2 dmu = <f,Af>
                   <= lambda^(-1) ||Af||_2^2.
```

Since `N(f) <= int|grad f|^2`, this proves `R_mu <= C_P(mu)`.

### Proof of the lower bound without an eigenfunction

Fix `epsilon>0`.  By the definition of the bottom of the spectrum, the
spectral projection of `A` on `[lambda,lambda+epsilon]` is nonzero.  Choose a
mean-zero unit vector `f_epsilon` in its range.  It belongs to `Dom(A)` and

```
int |grad f_epsilon|^2 dmu >= lambda,
||A f_epsilon||_2 <= lambda+epsilon.           (3)
```

For every coordinate function `x_j`, the form identity gives

```
(m_f)_j = q(f,x_j) = <Af,x_j>.
```

Consequently

```
|m_f| = sup_{|v|=1} |<Af,<v,x>>|
      <= ||Af||_2 sup_{|v|=1} ||<v,x>||_2
      = ||Af||_2,                              (4)
```

where the last equality is exactly isotropy.  Equations (3)--(4) imply

```
N(f_epsilon) >= lambda-(lambda+epsilon)^2.
```

Substitution in (1), followed by `epsilon downarrow 0`, yields

```
R_mu >= [lambda-lambda^2]/lambda^2
      = 1/lambda-1 = C_P(mu)-1.
```

This proof explicitly handles a continuous spectral edge.  In particular, it
does not assume that the Poincare infimum is attained.

## 2. Exact small-eigenfunction consequences

If the bottom is attained and `f` is normalized by

```
int f dmu=0,  int f^2 dmu=1,  Af=lambda f,
```

then (4) sharpens to

```
m_f = lambda int f(x)x dmu(x),
|m_f| <= lambda,                               (5)
```

and hence

```
int |grad f-m_f|^2 dmu >= lambda(1-lambda).    (6)
```

If `dmu=Z^(-1)e^(-V)dx`, with `V in C^2` convex and enough regularity to use
the weighted Bochner identity, then

```
||Af||_2^2
 = int ||Hess f||_HS^2 dmu
   + int Hess V(grad f,grad f) dmu
 = lambda^2.                                  (7)
```

Thus the normalized vector field `grad f/sqrt(lambda)` has centered `L^2`
mass at least `1-lambda`, while its full deformation-plus-curvature energy is
`lambda`.  Combining (6) and (7) with any universal estimate of the form

```
int |grad g-int grad g|^2
 <= C int [||Hess g||_HS^2+Hess V(grad g,grad g)]
```

would indeed prove `lambda >= 1/(C+1)`.  But (2) shows that the best possible
`C` in precisely this estimate is between `C_P-1` and `C_P`.  Calling it a
curl-free rigidity estimate does not make it weaker than KLS.

The Bochner identity is not needed to prove (2); it only identifies the
operator energy in (1) with the proposed geometric functional.  Therefore
weakening smoothness cannot remove the circularity.

## 3. Convex bodies and the Reilly boundary term

Let `K` be any bounded convex body with nonempty interior and let `mu_K` be
normalized Lebesgue measure on `K`.  The closed Neumann form

```
q_K(f,g)=|K|^(-1) int_K <grad f,grad g> dx,
Dom(q_K)=H^1(K)
```

defines the Neumann Laplacian `A_K`.  Convex bodies are Lipschitz domains, the
embedding `H^1(K) -> L^2(K)` is compact, and the first nonzero eigenvalue is
attained.  The proof of (2) applies verbatim to `A_K`: it uses only the form
identity with the coordinate functions.  In particular it already covers
polytopes and requires no boundary smoothing.

When `partial K` is `C^2`, is convex, and `f` is a smooth Neumann test
function, the Reilly identity (with the convention `II >= 0` on a convex
boundary) reads

```
||A_K f||_2^2
 = |K|^(-1) int_K ||Hess f||_HS^2 dx
   + |K|^(-1) int_{partial K}
       II(grad_T f,grad_T f) dH^(n-1).         (8)
```

For a smooth log-concave weight on a smooth convex domain, one adds

```
int_K Hess V(grad f,grad f) dmu
```

to the right side.  Therefore the best constant in the corresponding
Hessian-plus-curvature-plus-boundary rigidity estimate again lies in
`[C_P-1,C_P]`.

For a polytope, every open facet has `II=0`; the operator formulation remains
valid and is the rigorous nonsmooth version.  In particular, a proof that
requires a pointwise positive boundary-curvature contribution cannot cover
cubes or simplices uniformly.  Formula (8) may be recovered on smooth
approximants, but no such approximation is needed for the equivalence (2).

## 4. Model checks

### Gaussian

For the standard Gaussian, `C_P=1`, so (2) only gives `0 <= R_mu <= 1`.
The first eigenspace consists of affine functions and is annihilated by the
centering in `N(f)`.  A second-chaos Hermite polynomial gives ratio `1/2` in
(1), and the Hermite decomposition shows `R_mu=1/2`.  This is consistent with
(2) and explains the necessary subtractive `-1` in its lower bound.

### Isotropic interval and cube

On `[-sqrt(3),sqrt(3)]`, the first Neumann eigenfunction is

```
f(x)=sqrt(2) sin(pi x/(2sqrt(3))),
lambda=pi^2/12,
C_P=12/pi^2.
```

Its mean derivative has squared norm `2/3`, and its ratio in (1) is

```
12(pi^2-8)/pi^4,
```

a positive universal number.  Tensoring gives the same check on the
isotropic cube in every dimension.  On the cube all classical boundary
curvature terms vanish on the facets, so no hidden positive-curvature gain is
present.

### Exponential measures and non-attainment

For the centered one-sided exponential law

```
dmu(x)=exp(-(x+1)) 1_{x>=-1} dx,
```

the variance is one and the sharp Poincare constant is `4`.  The spectral gap
`lambda=1/4` is the bottom of continuous spectrum and is not an `L^2`
eigenvalue.  The spectral-band proof above still gives

```
3 <= R_mu <= 4.
```

The same example can be tensorized in arbitrary dimension without changing
`C_P`.  It is also strongly non-symmetric.  Thus neither tensor products,
asymmetry, unbounded support, nor lack of an eigenfunction creates a loophole
in the equivalence.

For the isotropic symmetric Laplace law, `C_P=2` and (2) gives
`1 <= R_mu <= 2`; its gap is likewise a continuous spectral edge.

## 5. Two further invalid reductions based on eigenfunctions

These examples do not use a hypothetical KLS counterexample.

1. **The eigenfunction pushforward need not be log-concave.**  If `X` is
   uniform on `[-a,a]` and `f(X)=sin(pi X/(2a))`, then `f(X)` has density

   ```
   rho(y)=1/[pi sqrt(1-y^2)],  -1<y<1.
   ```

   Since

   ```
   (log rho)''(y)=(1+y^2)/(1-y^2)^2 > 0,
   ```

   the density is strictly log-convex, not log-concave.  One-dimensional
   log-concave isoperimetry cannot be applied to `f_#mu`.

2. **Even a first nodal domain need not be convex.**  Put
   `a=sqrt(3)` and, on the isotropic square `[-a,a]^2`, let

   ```
   F(x,y)=2 sin(pi x/(2a))+sin(pi y/(2a)).
   ```

   Both summands lie in the first Neumann eigenspace, so `F` is a first
   eigenfunction.  The points

   ```
   p=(a/3,-a),  q=(0,0)
   ```

   satisfy `F(p)=F(q)=0`, whereas

   ```
   F((p+q)/2)=(sqrt(6)-2sqrt(2))/2 < 0.
   ```

   Hence the closed positive nodal set is not convex.  Small inward
   perturbations give the same conclusion for the open nodal domain.  Taking
   a product with isotropic intervals yields the counterexample in every
   dimension `n>=2`.

## 6. Audit conclusion

The exact information supplied by a small eigenvalue is (5)--(7).  Turning
that information into a universal lower bound by controlling the centered
gradient with its Bochner/Reilly energy requires a constant `R_mu` satisfying
(2), and hence requires the KLS bound itself.  Nodal convexity and
log-concavity of the eigenfunction distribution are false even for an
interval or a square.  This spectral-rigidity family is therefore blocked
unless one introduces an additional functional or geometric mechanism not
reducible to (1).
