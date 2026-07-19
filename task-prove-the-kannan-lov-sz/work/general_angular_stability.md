# General angular stability from Brenier splitting and Gaussian halfspace flux

## 1. Statement

Fix `delta in (0,1/2)`. Let `pi` be a `t`-strongly log-concave probability
on its affine hull, let `S` be Borel with

\[
 g=\pi(S)\in[\delta,1-\delta],
\]

and write, with `m=E_pi X`,

\[
 v=E[(1_S-g)(X-m)],\quad u={v\over|v|},\quad
 P=I-uu^T,
\]

\[
 D=E[(1_S-g)(X-m)(X-m)^T],\qquad
 \varepsilon=1-{\sqrt t|v|\over I(g)}.              \tag{1.1}
\]

For `v=0` the near-equality regime is absent. There is a modulus
`Omega_delta:[0,1]->[0,infinity)` tending to zero at zero such that

\[
\boxed{
 ||PD||_{HS}^2\le {C_\delta\over t^2}
       \Omega_\delta(\varepsilon),}                  \tag{1.2}
\]

and one may take

\[
 \rho_\delta(r)=C_\delta\left\{\sqrt r+
       \big(r\sqrt{\log(e/r)}\big)^{1/3}\right\},
\qquad
 \Omega_\delta(r)=r+\sqrt{\rho_\delta(r)}.           \tag{1.3}
\]

Thus, up to a logarithmic factor, `Omega_delta(r)=O_delta(r^(1/6))`.
No ambient-dimensional constant occurs.

## 2. Standardization and a Gaussian halfspace in the source

It is enough to prove the assertion for `t=1`, `m=0`; restore scale at the
end.  All Gaussian spaces, matrices, and Hilbert--Schmidt norms in the proof
are taken on the parallel linear space of the affine hull of `pi`.  Let `G`
be standard Gaussian there and let `T=grad Phi` be the centered
Brenier contraction taking `G` to `X~pi`. Its a.e. Jacobian is a symmetric
positive contraction.

The clean-room spatial stability theorem gives two facts. First, if `H` is
the active target halfspace of mass `g`, then

\[
                         \pi(S\mathbin\triangle H)
 \le C_\delta\sqrt\varepsilon.                       \tag{2.1}
\]

Second, for the active scalar coordinate,

\[
 E|\langle u,T(G)\rangle-\langle u,G\rangle|^2
 \le \zeta_\delta(\varepsilon),\qquad
 \zeta_\delta(r)=C_\delta r\sqrt{\log(e/r)}.         \tag{2.2}
\]

Let `a` be the threshold defining `H` in the target coordinate and let
`b=Phi^{-1}(1-g)`. The elementary threshold-smoothing lemma says that if
`Y,Z` have `E|Y-Z|^2<=zeta`, `Z` is standard Gaussian, and the two upper
threshold events have the same central mass `g`, then

\[
 P(1_{Y\ge a}\ne1_{Z\ge b})\le C_\delta\zeta^{1/3}. \tag{2.3}
\]

Indeed, put `E={Y>=a}`, `F_a={Z>=a}`, and `F={Z>=b}`.  Chebyshev and the
bounded Gaussian density give

\[
 P(E\mathbin\triangle F_a)\le \zeta/h^2+2\varphi(0)h .
\]

The Gaussian threshold events `F_a,F` are nested and `P(E)=P(F)`, hence
`P(F_a triangle F)<=P(E triangle F_a)`.  The triangle inequality and
`h=zeta^(1/3)` prove (2.3), without first estimating `|a-b|`.

Put

\[
 A=T^{-1}(S),\qquad B=\{\langle G,u\rangle\ge b\}.
\]

Both have Gaussian mass `g`. Equations (2.1)--(2.3) give

\[
\boxed{\gamma_n(A\mathbin\triangle B)\le
       \rho:=\rho_\delta(\varepsilon).}              \tag{2.4}
\]

## 3. Dimension-free perturbation of the first two correlated moments

Let

\[
 \begin{aligned}
 v_A&=E[(1_A-g)X],&D_A&=E[(1_A-g)XX^T],\\
 v_B&=E[(1_B-g)X],&D_B&=E[(1_B-g)XX^T],
 \end{aligned}                                      \tag{3.1}
\]

where `X=T(G)` is centered. Since `Cov(X)<=I`, Cauchy--Schwarz in every
direction gives

\[
\boxed{|v_A-v_B|\le\sqrt\rho.}                       \tag{3.2}
\]

The quadratic perturbation is also dimension free. For every symmetric
matrix `M` with `||M||HS=1`, the Poincare inequality of the one-strongly
log-concave target gives

\[
 Var(X^TMX)\le E|2MX|^2\le4.                          \tag{3.3}
\]

Because `E(1_A-1_B)=0`, (3.3) and Cauchy--Schwarz imply

\[
\boxed{||D_A-D_B||_{HS}\le2\sqrt\rho.}               \tag{3.4}
\]

The same calculation with `M=theta theta^T` gives

\[
\boxed{||D_B||_{op}\le1.}                            \tag{3.5}
\]

All three estimates remain valid for hard support by the closed Poincare
form and approximation; no trace or `sqrt(n)` enters.

## 4. Apply the exact halfspace-pullback theorem

Let `u_B=v_B/|v_B|` and `P_B=I-u_Bu_B^T`. For sufficiently small
`epsilon`, (3.2) and centrality ensure `v_B!=0`, and

\[
 |u-u_B|\le C_\delta\sqrt\rho.                       \tag{4.1}
\]

The halfspace-pullback contraction theorem applies to the fixed Gaussian
halfspace `B` and the arbitrary nonlinear Brenier contraction `T`:

\[
 ||P_BD_B||_{HS}^2
 \le C_\delta\{I(g)-|v_B|\}.                         \tag{4.2}
\]

Moreover,

\[
 I(g)-|v_B|
 \le I(g)-|v_A|+|v_A-v_B|
 \le C_\delta\varepsilon+\sqrt\rho.                 \tag{4.3}
\]

Using (3.4), (3.5), and (4.1),

\[
 \begin{aligned}
 ||P D_A||_{HS}
 &\le ||D_A-D_B||_{HS}
   +||(P-P_B)D_B||_{HS}+||P_BD_B||_{HS}\\
 &\le C_\delta\left\{\sqrt\rho+
          (\varepsilon+\sqrt\rho)^{1/2}\right\}.
 \end{aligned}                                      \tag{4.4}
\]

Here `P-P_B` has rank at most two, so
`||(P-P_B)D_B||HS<=||P-P_B||HS||D_B||op`; this is the point at which an
incorrect trace estimate would introduce dimension. Squaring (4.4) proves

\[
 ||PD_A||_{HS}^2\le C_\delta
       \{\varepsilon+\sqrt\rho\}.                    \tag{4.5}
\]

For large `epsilon`, enlarge the constant and use the universal quadratic
Poincare bound, so (4.5) holds on all of `[0,1]` with the modulus in (1.3).

Finally, under the standardization `X_tilde=sqrt(t)(X-m)`, one has
`D_tilde=tD`, while the projection is unchanged. Thus (4.5) becomes (1.2).

## 5. Consequences and limitation

At a central near-equality posterior, `|v| asymp_delta t^(-1/2)`, so (1.2)
gives the fixed-time tilt derivative estimate

\[
 {||PD||_{HS}^2\over|v|^2}
 \le {C_\delta\over t}\Omega_\delta(\varepsilon).   \tag{5.1}
\]

This supplies the previously missing pointwise transverse control for an
arbitrary set and arbitrary strongly log-concave posterior. It composes only
already proved statements and is compatible with hard support and
lower-dimensional affine hulls.

The modulus is too weak by itself to globalize phase directions across a
long localization window: selecting extremely small `epsilon` by a profile
pigeonhole can enlarge the time/overlap scale faster than (1.3) improves.
Thus (1.2) is a genuine new lemma, not a completed KLS proof. Its next use
must combine fixed-scale phase expansion, temporal defect dissipation, or
the calibrated-ray endpoint constraints without reintroducing that scale
mismatch.
