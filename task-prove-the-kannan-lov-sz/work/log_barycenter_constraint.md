# Logarithmic barycenters of small events and coherent ray packets

## 1. A sharp small-event barycenter bound

Let `mu` be isotropic and log-concave on `R^n`.  If `E` is a Borel event of
mass `epsilon in (0,1)`, then

\[
 \boxed{\quad
 \left|\mathbb E[X\mid E]\right|
 \le C\left(1+\log {1\over\epsilon}\right).
 \quad}                                                           \tag{1}
\]

Indeed, let `b=E[X|E]` and, when `b!=0`, put `theta=b/|b|`.  The marginal
`Y=<theta,X>` is a centered variance-one log-concave random variable.  Its
one-dimensional exponential-tail estimate is

\[
 \mathbb P\{Y\ge t\}\le C_0e^{-c_0t}\qquad(t\ge0).                \tag{2}
\]

Among all events of mass `epsilon`, the conditional mean of `Y` is maximized
by an upper level set.  Integrating (2), or integrating the decreasing
quantile of `Y`, gives

\[
 \mathbb E[Y\mid E]
 \le C\left(1+\log {1\over\epsilon}\right).                       \tag{3}
\]

The left side is `|b|`, proving (1).  This improves the covariance-only
bound `|b|<=epsilon^{-1/2}` whenever `epsilon` is small.

## 2. Application to balanced transport rays

Let a probability be disintegrated over oriented affine rays:

\[
 X=z_y+T N_y,\qquad |N_y|=1,qquad d\mu=d\mu_y(T)d\eta(y).          \tag{4}
\]

Assume that every ray has the same sign proportions

\[
 \mu_y(T>0)=p,qquad \mu_y(T<0)=q=1-p,                             \tag{5}
\]

where `p,q>=delta`.  Put

\[
 d_y=\mathbb E_y[T\mid T>0]-\mathbb E_y[T\mid T<0]>0.             \tag{6}
\]

For a quotient event `Omega_i`, define

\[
 B_i=\{y\in\Omega_i,T>0\},
 \qquad C_i=\{y\in\Omega_i,T<0\}.                                \tag{7}
\]

Because the sign proportions in (5) are independent of `y`, conditioning
on `B_i` and `C_i` induces the same normalized quotient law
`eta(.|Omega_i)`.  The base points cancel exactly, and hence

\[
 b_{B_i}-b_{C_i}
 =\mathbb E_{\eta(\cdot|\Omega_i)}[d_yN_y].                       \tag{8}
\]

The two events have masses `p eta(Omega_i)` and `q eta(Omega_i)`.  Applying
(1) twice gives

\[
 \boxed{\quad
 \left|\mathbb E_{\eta(\cdot|\Omega_i)}[d_yN_y]\right|
 \le C_\delta\left(1+\log {1\over\eta(\Omega_i)}\right).
 \quad}                                                           \tag{9}
\]

In particular, if `d_y>=s` on `Omega_i` and all its directions lie in the
unit Euclidean cap about a unit vector `u_i`, then
`<N_y,u_i>>=1/2`, so (9) yields

\[
 s\le C_\delta\left(1+\log {1\over\eta(\Omega_i)}\right).        \tag{10}
\]

Thus a bad scale `s` cannot be carried by a bounded number of coherent
direction packets of polynomially small mass.  Its only remaining escape is
an extremely diffuse direction law, with every coherent cap having mass
exponentially small in `s`.  This is exactly the branch in which a global
parallel-versus-concurrent (or focal-turning) inverse theorem is still
needed; (10) alone does not control a nearly uniform direction law in an
ambient dimension much larger than `s^2`.

