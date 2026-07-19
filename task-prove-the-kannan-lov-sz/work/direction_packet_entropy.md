# Exponential smallness of coherent long-ray packets

This note turns the logarithmic barycenter bound for small events into a
dimension-free restriction on every coherent packet of long balanced rays.

## 1. Barycenters of small events

**Lemma 1.**  Let `X` be isotropic and log-concave in `R^n`.  If
`E` is Borel and `P(E)=epsilon in (0,1/2]`, then

\[
 \left|\mathbb E[X\mid E]\right|
 \le C_b\bigl(1+\log(1/\epsilon)\bigr).                  \tag{1}
\]

**Proof.**  Put `b=E[X|E]`; the assertion is trivial if `b=0`.
For `theta=b/|b|`, the variable `Y=<X,theta>` is one-dimensional,
centered, variance one, and log-concave.  The standard one-dimensional
log-concave tail estimate gives universal constants `C_0,c_0>0` with

\[
 \mathbb P(Y\ge t)\le C_0e^{-c_0t}\qquad(t\ge0).          \tag{2}
\]

For every nonnegative random variable `Z` and every event of mass
`epsilon`, layer cake gives

\[
 \mathbb E[Z{\bf1}_E]
 \le\int_0^\infty\min\{\epsilon,\mathbb P(Z\ge t)\}\,dt. \tag{3}
\]

Apply (3) to `Y_+` and split the integral at
`t_0=c_0^{-1}\log(C_0/epsilon)`.  Equations (2)--(3) give

\[
 \epsilon|b|=\mathbb E[Y{\bf1}_E]
 \le\mathbb E[Y_+{\bf1}_E]
 \le \epsilon t_0+C_0c_0^{-1}e^{-c_0t_0}.
\]

This is (1), after changing the universal constant.  QED.

## 2. Application to balanced rays

Let

\[
 \mu=\int \nu_y\,d\eta(y),\qquad X=z_y+T N_y           \tag{4}
\]

be a nonbranching ray disintegration of an isotropic log-concave law, with
`|N_y|=1`.  Fix a quotient set `Omega`.  Assume that for every
`y in Omega` there are measurable lower and upper conditional tails
`L_y,U_y` such that

\[
 \nu_y(L_y)=\nu_y(U_y)=\kappa,
 \qquad
 \mathbb E[T\mid U_y]-\mathbb E[T\mid L_y]\ge d_0s,     \tag{5}
\]

where `kappa,d_0>0` are fixed universal constants.  Quantile tails make the
choice measurable; the constant-mass two-sided excursions of a true T3
extremizer provide (5).

For a measurable packet `P subset Omega`, put `epsilon=eta(P)` and

\[
 B_P=\{z_y+tN_y:y\in P,\ t\in U_y\},\qquad
 C_P=\{z_y+tN_y:y\in P,\ t\in L_y\}.                    \tag{6}
\]

Then `mu(B_P)=mu(C_P)=kappa epsilon`, and exact cancellation of the
ray bases gives

\[
 b_{B_P}-b_{C_P}
 ={1\over\epsilon}\int_P d_yN_y\,d\eta(y),
 \qquad d_y:=E[T|U_y]-E[T|L_y]\ge d_0s.                 \tag{7}
\]

**Proposition 2 (coherent packet bound).**  If there is a unit vector `v`
such that

\[
 \langle N_y,v\rangle\ge\gamma>0
 \quad\text{for eta-a.e. }y\in P,                       \tag{8}
\]

then

\[
 \boxed{\quad
  \epsilon\le {1\over\kappa}
       \exp\!\left(1-c\gamma d_0s\right),
 \quad}                                                   \tag{9}
\]

with a universal `c>0`.

**Proof.**  Taking the scalar product of (7) with `v` gives

\[
 |b_{B_P}-b_{C_P}|\ge\gamma d_0s.                        \tag{10}
\]

Lemma 1 applied to the two events in (6) gives

\[
 |b_{B_P}-b_{C_P}|
 \le2C_b\bigl(1+\log(1/(\kappa\epsilon))\bigr).          \tag{11}
\]

Rearranging (10)--(11) proves (9).  QED.

For example, every chordal unit cap
`{u:|u-v|<=1}` satisfies `<u,v>>=1/2`.  Thus, on a
positive-mass family of rays of scale `s`, every such direction cap has
quotient mass at most `C exp(-cs)`.  In particular the direction law has
unit-scale covering/entropy number at least

\[
 N_{\rm dir}\ge c\,e^{cs}.                               \tag{12}
\]

This improves the covariance-only lower bound `N_dir>=cs^2` to an
exponential one.

## 3. What remains

Proposition 2 kills every finite or polynomial-size family of coherent long
cylinders.  The only possible escape is an exponentially diffuse direction
law.  A concurrent radial family has exactly this qualitative behavior in
high dimension, so diffuseness alone is not a contradiction.  To finish the
inverse step one would need to combine (12) with normal congruence:

* smooth direction change is constrained by
  `s^2||S_y||_HS^2<=C`;
* discontinuous direction change must occur through focal/medial junctions;
* global log-concavity supplies bridges between the exponentially many caps.

The existing indicator midpoint count does not perform this last conversion:
it is null-set unstable and gives no mass to actual reflected
representations.  Thus (9) is a genuine new constraint, but an additional
weighted bridge or full-BV turning estimate is still required to identify the
diffuse branch with concurrence.
