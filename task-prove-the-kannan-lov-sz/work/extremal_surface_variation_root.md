# Formal second variation of a signed-distance T3 extremizer

This note records a smooth-chart calculation which is being independently
audited.  Its purpose is to extract information special to a *global
first-moment extremizer*, beyond balance of its transport rays.

## 1. Shape derivatives of signed distance

Let `Sigma={f=0}` be a smooth oriented hypersurface and suppose

\[
 F(z,t)=z+tN(z)
\]

is a nonsingular normal chart on the part under consideration.  Let `S=DN`
with the convention `D_zF=I+tS`.  Deform `Sigma` with normal speed `h` and
let `f_epsilon` be the signed distance to the deformed surface, with the same
orientation.  At a fixed point `F(z,t)`, write derivatives at zero as dots.
Then

\[
 \dot f(z,t)=-h(z),                                      \tag{1}
\]

and

\[
 \boxed{\quad
 \ddot f(z,t)=-\int_0^t
 \left|(I+sS_z)^{-1}\nabla_\Sigma h(z)\right|^2ds.
 \quad}                                                  \tag{2}
\]

Indeed, differentiating `|grad f_epsilon|^2=1` once shows that `dot f` is
constant on each original normal ray, and its boundary value is `-h`.
Differentiating a second time gives

\[
 \partial_t\ddot f=-|\nabla\dot f|^2
 =-|(I+tS)^{-1}\nabla_\Sigma h|^2.
\]

The twice-differentiated boundary condition gives `ddot f(z,0)=0`, proving
(2).  For a tilted plane, (2) reduces to the exact expansion
`ddot f=-t|grad h|^2`.

## 2. Balance kills every first variation

Disintegrate `mu` in the chart as

\[
 d\mu=q_y(t)dt\,d\eta(y),
\]

and assume

\[
 \mu_y(t>0)=p,\qquad\mu_y(t<0)=q=1-p                    \tag{3}
\]

for almost every ray.  Put `s_0=p-q`.  If `bar h=int h deta`, then

\[
 {d\over d\epsilon}\bigg|_0
 \left(f_\epsilon-\int f_\epsilon d\mu\right)
 =-h+\bar h.                                            \tag{4}
\]

Since `E_y sign(T)=s_0` on every ray, (4) gives

\[
 {d\over d\epsilon}\bigg|_0
 \int|f_\epsilon-\mu f_\epsilon|d\mu=0                 \tag{5}
\]

for every compactly supported smooth `h`.  Thus the constant sign
proportions are exactly the Euler--Lagrange equation for normal displacement.

## 3. The extremal surface-stability inequality

Assume now that `f` is a global maximizer of the centered `L^1` objective
over all 1-Lipschitz functions, and that the deformation remains within a
nonsingular chart.  The distributional second derivative of absolute value
gives

\[
 \begin{split}
 J''(0)
 &=\mathbb E[(\operatorname {sign}T-s_0)\ddot f]\\
 &\quad+2\int q_y(0)(h(y)-\bar h)^2d\eta(y).             \tag{6}
 \end{split}
\]

For `t>0`, `sign(t)-s_0=2q` and (2) is nonpositive.  For `t<0`,
`sign(t)-s_0=-2p` and (2) is nonnegative.  Since maximality gives
`J''(0)<=0`, (6) yields

\[
 \boxed{\begin{split}
 &\int q_y(0)(h-\bar h)^2d\eta(y)\\
 &\le q\int d\eta(y)\int_{t>0}q_y(t)
       \int_0^t |(I+sS_y)^{-1}\nabla h|^2ds\,dt\\
 &\quad+p\int d\eta(y)\int_{t<0}q_y(t)
       \int_t^0 |(I+sS_y)^{-1}\nabla h|^2ds\,dt.
                                                               \tag{7}
 \end{split}}\]

Equivalently, after Tonelli, the right side is a positive weighted Dirichlet
form on the zero surface, with the positive and negative conditional tail
functions as coefficients.

For a flat product chart, (7) is a genuine quotient Poincare inequality with
boundary weight `q_y(0)` and ray first-moment weight on the energy.  This is
consistent with tensorization: if the quotient had a worse scale, the flat
coordinate could not be the global T3 extremizer.

## 4. Normal-coordinate test

For `h_a(z)=<a,N(z)>`,

\[
 \nabla_\Sigma h_a=S_zP_{T_z\Sigma}a.
\]

Summing (7) over an orthonormal ambient basis gives on the left the weighted
normal variance

\[
 \int q_y(0)|N_y-m_0|^2d\eta(y),
 \qquad
 m_0={\int q_y(0)N_yd\eta(y)\over\int q_y(0)d\eta(y)},   \tag{8}
\]

up to the distinction between the unweighted centering `bar h` in (7) and
the boundary-weighted minimizer in (8), which only lowers the expression.
The summed integrand on the right is

\[
 \operatorname {Tr}\left[S_y^2(I+sS_y)^{-2}\right].      \tag{9}
\]

The polynomial-Jacobian lemma gives

\[
 \sigma_y^2|S_y|_{HS}^2\le C_\delta                    \tag{10}
\]

on every smooth ray whose two sign masses are bounded below.  Equations
(7)--(10) therefore couple dispersion of the normals to the *full* shape
curvature at the correct ray scale.  They do not yet give normal alignment:
the Gaussian fan carries the missing variation at singular focal junctions,
where (2) is not a complete description of the deformation.

## 5. Remaining rigor and the target use

A final use of (7) requires all of the following:

1. justification of signed-distance shape derivatives up to the finite ray
   endpoints and passage through the centering operation;
2. approximation of a nonsmooth global extremizer without losing the exact
   maximality inequality;
3. a second-variation term on the completed normal cycle which charges
   polyhedral, medial, and focal junctions; and
4. a proof that the combined smooth-plus-singular form forces either a
   parallel factor or an approximately concurrent family.

No quotient Poincare inequality may be inserted in item 4: on a flat product
chart that would simply reintroduce the lower-dimensional KLS constant.

