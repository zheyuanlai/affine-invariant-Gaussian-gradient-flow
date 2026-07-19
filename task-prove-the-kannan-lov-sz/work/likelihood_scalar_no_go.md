# A scalar likelihood weight cannot repair normalized-localization transfer

This note isolates a limitation of endpoint functionals which multiply the
entire anisotropic posterior perimeter by a nonnegative function of the two
global statistics

```
g_t=mu_t(S),
r_t={v_t^T A_t^{-1}v_t}/{[g_t(1-g_t)]},
v_t=Cov_t(1_S,X).
```

The limitation is face-specific: a harmless tail face can have arbitrarily
large posterior curvature weight while `(g_t,r_t)` converges to the state of
an unrelated transverse cut.  Consequently a useful endpoint functional has
to distinguish which boundary face carries the posterior uncertainty.  A
scalar function of `(g_t,r_t)` cannot do so.

## 1. What the projection ratio actually says

Let

```
h=1_S-g,
ell(x)=v^T A^{-1}(x-a).
```

Orthogonal projection in `L^2(mu)` gives the exact identity

```
||h-ell||_2^2=g(1-g)(1-r).                         (1)
```

Thus only a regime `r -> 1` makes the indicator asymptotically linear.  The
universal projection-gap lemma in `work/correlation_lemmas.md` gives

```
r <= 1-eta                                                   (2)
```

for a numerical `eta>0` for every log-concave posterior.  Hence the regime
of arbitrarily accurate linear prediction is actually absent.  A fixed
threshold `r>=r_0` gives only a fixed-error approximation in (1).  For the
uniform law on an interval cut at its median, `r=3/4`, so even the most basic
one-dimensional halfspace leaves one quarter of the Bernoulli variance in
the nonlinear residual.

## 2. Tail-decorated product cut

Let `eta_lambda` be the variance-one symmetric Laplace law on the line,

```
rho_lambda(s)=(lambda/2) exp(-lambda |s|),
lambda=sqrt(2).
```

Let `nu` be an isotropic log-concave probability on `R^d`, and let `B` be a
finite-perimeter set with `nu(B)=1/2`.  For `L>0`, put

```
D_L = ((-infinity,L] x B)
      union ((-infinity,-L] x B^c).                 (3)
```

Write `p=P_nu(B)` and `alpha_L=1-exp(-lambda L)`.  The mass and ordinary
perimeter are

```
(eta_lambda tensor nu)(D_L)=1/2,
P_{eta_lambda tensor nu}(D_L)=rho_lambda(L)+alpha_L p.        (4)
```

Under covariance-normalized localization, product measures stay products.
Denote the two posterior factors by `eta_{L,t}` and `nu_t`, their accumulated
precisions by `q_t` and `Q_t^nu`, and put

```
b_t=nu_t(B),
s_t={Cov_{nu_t}(1_B,Y)^T (A_t^nu)^{-1}
                         Cov_{nu_t}(1_B,Y)}
       /[b_t(1-b_t)].                                      (5)
```

The likelihood law at a signal `(x,z)` also factors as
`P^(x,z)=P_eta^x tensor P_nu^z`.

For a bounded continuous function
`Phi:(0,1)x[0,1] -> [0,infinity)`, define the scalar-weighted endpoint
perimeter

```
R_{T,Phi}(D)
 = E[ Phi(g_T,r_T) P_{Q_T^{-1},mu_T}(D) ].            (6)
```

The two horizontal faces of `D_L` will be denoted by `H_L^+` and `H_L^-`.

## 3. Exact horizontal-face asymptotic

### Proposition

Fix `T>0`.  Suppose in addition that `Phi` vanishes when its first argument
lies outside a compact subinterval of `(0,1)`.  Then

```
lim_{L->infinity}
  R_{T,Phi}^{horizontal}(D_L)
    /[rho_lambda(L) sqrt(L)]
 = [lambda exp(T)]^{-1/2} E[Phi(b_T,s_T)].            (7)
```

The expectation on the right is for the normalized localization of the
transverse pair `(nu,B)` under the original path law.  For a fixed bounded
uniformly continuous `Phi` with the stated support, the convergence is
uniform over the choice of the transverse pair `(nu,B)`.  Indeed, every
error below is bounded by the first-factor probability in (9) and the
modulus of continuity of `Phi`; no transverse moment enters the estimate.

### Proof

Consider first a point `(L,z)` on `H_L^+`, so `z` is in the
measure-theoretic interior of `B`.  Under `P_eta^L`, Section 8 of
`work/normalized_localization.md` proves

```
q_T^{-1/2}/sqrt(L) -> [lambda exp(T)]^{-1/2}          (8)
```

in `L^1`; the lower-tail estimate there gives the required uniform
integrability.

The same posterior asymptotics show

```
eta_{L,T}((-L,L)) -> 1                               (9)
```

in probability and hence in `L^1`.  For completeness, the posterior mode is
`L(1-exp(-T))+o(L)`, its variance is `L/[lambda exp(T)]+o(L)`, and its
distance from either endpoint of `(-L,L)`, measured in posterior standard
deviations, tends to infinity.  The log-concave one-dimensional tail bound
then gives (9).

Let `Z,Y` be conditionally independent samples from the two posterior
factors.  From (3),

```
|1_{D_L}(Z,Y)-1_B(Y)| <= 1_{|Z|>=L}.                 (10)
```

Therefore the two centered indicators differ in `L^2` by `o(1)`.  Orthogonal
projection is a contraction in `L^2`.  Since the centered coordinate span
for the product is the orthogonal sum of the centered spans of its two
factors, (10) implies

```
g_T(D_L)-b_T -> 0,
v_T(D_L)^T A_T^{-1}v_T(D_L)
 - Cov_{nu_T}(1_B,Y)^T(A_T^nu)^{-1}Cov_{nu_T}(1_B,Y) -> 0      (11)
```

in probability.  On the compact `b_T`-range selected by `Phi`, division by
the Bernoulli variance is uniformly continuous, and hence

```
Phi(g_T(D_L),r_T(D_L))-Phi(b_T,s_T) -> 0             (12)
```

in probability.  Boundedness of `Phi`, (8), uniform integrability, and the
independence of the two posterior factors now give

```
E^(L,z)[Phi(g_T,r_T) q_T^{-1/2}]/sqrt(L)
 -> [lambda exp(T)]^{-1/2} E_nu^z[Phi(b_T,s_T)].     (13)
```

The positive horizontal face has normal `e_1`.  The posterior anisotropic
surface weight there is exactly `q_T^{-1/2}`.  The likelihood identity and
Fubini turn its contribution into

```
rho_lambda(L) int_B
 E^(L,z)[Phi(g_T,r_T)q_T^{-1/2}] dnu(z).             (14)
```

The negative face is identical, with `z` in `B^c` and signal `-L`.
Integrating (13) and its negative-face analogue uses the exact
change-of-measure identities

```
int_B E_nu^z[F]dnu(z)=E[b_T F],
int_{B^c} E_nu^z[F]dnu(z)=E[(1-b_T)F].               (15)
```

The two coefficients add to `E[Phi(b_T,s_T)]`.  Equations (13)--(15) prove
(7).

## 4. Quantified consequence for a hypothetical bad sequence

Let `(nu_j,B_j)` be any sequence of isotropic log-concave half-mass cuts
with `p_j=P_{nu_j}(B_j)->0`.  Fix `alpha` in `(0,1/2)` and choose `L_j` by

```
rho_lambda(L_j)/p_j=[log(1/p_j)]^{-alpha}.           (16)
```

Then `L_j~lambda^{-1}log(1/p_j)`, and (4) gives

```
P_{eta_lambda tensor nu_j}(D_{L_j})=(1+o(1))p_j.     (17)
```

Thus the decoration preserves asymptotic near-optimality.  If a family of
scalar weights satisfies

```
E[Phi(b_{j,T},s_{j,T})] >= kappa_T>0                 (18)
```

uniformly along the sequence, (7) and (16) yield

```
R_{T,Phi}^{horizontal}(D_{L_j})
 /P_{eta_lambda tensor nu_j}(D_{L_j})
 >= c_{T,lambda} kappa_T
       [log(1/p_j)]^{1/2-alpha} -> infinity.         (19)
```

Consequently there is no universal transfer inequality

```
R_{T,Phi}(D) <= C P_mu(D)                            (20)
```

for a scalar weight which remains nontrivial on the posterior states of a
putative bad core.  This is a proof-search no-go, not an assumption that a
bad sequence exists: in a contradiction proof of KLS, the tail-decorated
sequence is available whenever the sequence being excluded is available.

## 5. Balanced-mass and projection-ratio weights are covered

For a half-mass cut, the universal projection gap improves the posterior
mass estimate to

```
E[b_T(1-b_T)] >= exp(-(1-eta)T)/4.                  (21)
```

Since `min(b,1-b)>=b(1-b)`, put

```
m_T=exp(-(1-eta)T)/4.
```

Choose a continuous cutoff `chi_T:[0,1]->[0,1]` which is one on
`[m_T/2,1-m_T/2]` and zero outside `[m_T/4,1-m_T/4]`.  From

```
E min(b_T,1-b_T)>=m_T
```

one obtains

```
P{m_T/2<=b_T<=1-m_T/2}>=m_T,
E chi_T(b_T)>=m_T.                                  (22)
```

(The weaker positive universal lower bound is all that is needed.)  Hence
`Phi(g,r)=chi_T(g)` satisfies (18).  Splitting this cutoff into any
nonnegative low- and high-projection pieces,

```
Phi_low(g,r)+Phi_high(g,r)=chi_T(g),                 (23)
```

does not help: the sum of the two horizontal contributions still has the
divergence (19), and at least one branch has at least half of it.

This rules out the proposed repair in which one retains the full
anisotropic perimeter and merely weights or stops it according to the
global pair `(g_t,r_t)`.  The obstruction is that `(g_t,r_t)` describes the
transverse core while the amplified boundary normal belongs to the
independent tail coordinate.  A viable likelihood-aware functional must be
face-specific--for example a signed calibration or an influence/flux whose
normal component vanishes on posterior-irrelevant faces--rather than a
nonnegative scalar multiple of the full endpoint perimeter.
