# Spectral-window Bochner ledger at a non-attained gap

## Status and purpose

This note replaces every use of a first Poincare eigenfunction in the
near-linear Bochner calculations by a bottom spectral window.  It is an
unconditional operator lemma.  In particular it applies when the bottom of
the nonconstant spectrum is continuous, as it is for the centered
one-sided exponential and the symmetric Laplace laws.

The note does **not** prove KLS.  Its purpose is to remove spectral
attainment as a loophole in the existing rigidity audits.

## 1. Closed-form setting

Let \(\mu\) be a centered isotropic log-concave probability on
\(\mathbb R^n\), initially with full-dimensional support.  Let

\[
 \mathcal E(u,v)=\int\langle\nabla u,\nabla v\rangle\,d\mu
\]

be the closed Dirichlet form, let \(A\geq0\) be its self-adjoint generator,
and put

\[
 \lambda=\inf\sigma(A|_{\mathbf 1^\perp})=C_P(\mu)^{-1}.
\]

For every \(\varepsilon>0\), the spectral projection
\(\mathbf 1_{[\lambda,\lambda+\varepsilon]}(A)\) on
\(\mathbf1^\perp\) is nonzero.  Choose a unit vector \(f\) in its range and
write

\[
 q=\langle Af,f\rangle=\lambda+\alpha,
 \qquad z=(A-\lambda)f,
 \qquad 0\leq\alpha\leq\varepsilon .                 \tag{1.1}
\]

The spectral theorem gives the sharper window-error estimate

\[
 \boxed{\ \|z\|_2^2\leq\varepsilon\alpha\leq\varepsilon^2.\ } \tag{1.2}
\]

Indeed, on the window the spectral variable \(t=s-\lambda\) belongs to
\([0,\varepsilon]\), so \(t^2\leq\varepsilon t\).

## 2. Affine projection and mean-gradient error

Set

\[
 a=\int x f(x)\,d\mu(x),\qquad r=f-a\cdot x .          \tag{2.1}
\]

Isotropy says that \(a\cdot x\) is exactly the \(L^2(\mu)\)-orthogonal
projection of \(f\) onto the linear functions.  Thus

\[
 \int r\,d\mu=0,\qquad \int xr\,d\mu=0,
 \qquad \|r\|_2^2=1-|a|^2.                            \tag{2.2}
\]

For every coordinate function, the weak form identity yields

\[
 \int\partial_i f\,d\mu=\langle Af,x_i\rangle.
\]

Consequently, with

\[
 d=\int xz(x)\,d\mu(x),
\]

one has

\[
 \boxed{\ \int\nabla f\,d\mu=\lambda a+d,\qquad |d|\leq\|z\|_2
 \leq\sqrt{\varepsilon\alpha}\leq\varepsilon.\ }    \tag{2.3}
\]

The vector estimate follows from isotropy: for every centered scalar
\(g\),

\[
 \left|\int xg\,d\mu\right|
 =\sup_{|u|=1}\left|\int (u\cdot x)g\,d\mu\right|
 \leq\|g\|_2.
\]

In particular,

\[
 \operatorname {Var}_\mu(\nabla f)
 =\lambda+\alpha-|\lambda a+d|^2.                    \tag{2.4}
\]

## 3. The corrected nonnegative Bochner ledger

Assume first that \(d\mu=Z^{-1}e^{-V}1_\Omega dx\), where \(V\) and the
convex domain \(\Omega\) are smooth enough for the weighted
Bochner--Reilly identity.  Denote

\[
 B=\int\|D^2f\|_{\rm HS}^2\,d\mu
\]

and by \(\mathcal R\geq0\) the sum of the potential-curvature and convex
boundary-curvature terms.  Then

\[
 B+\mathcal R=\|Af\|_2^2
 =\lambda^2+2\lambda\alpha+\|z\|_2^2.                \tag{3.1}
\]

Componentwise Poincare gives

\[
 B\geq\lambda\operatorname {Var}_\mu(\nabla f).       \tag{3.2}
\]

Substituting (2.4) into (3.1) gives the exact identity

\[
 \boxed{
 [B-\lambda\operatorname {Var}_\mu(\nabla f)]+\mathcal R
 =\lambda\alpha+\|z\|_2^2+\lambda|\lambda a+d|^2 .}
                                                               \tag{3.3}
\]

Every term on both sides of (3.3) is nonnegative.  This is the appropriate
replacement for the attained-eigenfunction ledger.  It also gives

\[
 \|D^2r\|_2=\|D^2f\|_2\leq\|Af\|_2\leq\lambda+\varepsilon.
                                                               \tag{3.4}
\]

The same statements hold on a smooth convex body with the Neumann
generator, with the Reilly boundary term included in \(\mathcal R\).  For a
nonsmooth log-concave law, (1.1)--(2.4) remain literal closed-form
identities.  Formula (3.3) is used only after the usual smooth convex
approximation, or equivalently in the weak Bochner formulation; no estimate
in this note discards a boundary term.

## 4. Passage to the spectral edge

Let \(\varepsilon_j\downarrow0\) and choose the above unit spectral-window
vectors.  Then

\[
 q_j\to\lambda,\qquad \|z_j\|_2\to0,
 \qquad \left|\int\nabla f_j-\lambda a_j\right|\to0,
 \qquad \|D^2r_j\|_2\leq\lambda+o(1).                \tag{4.1}
\]

Thus every contradiction argument based only on the exact first-eigenfunction
relations must also exclude these minimizing sequences.  The lack of an
attained eigenfunction cannot itself close the rigidity gap.

For example, the variance-one symmetric Laplace law has
\(C_P=2\), so \(\lambda=1/2\) is a non-attained continuous spectral edge.
The centered variance-one one-sided exponential has \(C_P=4\), so its edge
is \(1/4\).  Both are covered by (4.1).

## 5. Audit conclusion

The window error is genuinely \(O(\varepsilon)\), and all constants above
are dimension free.  The load-bearing unresolved statements in the
near-linear route (notably the residual mean-gradient estimate and the
nonlinear-mode exclusion) remain unresolved after this replacement; they
can no longer be bypassed by assuming spectral attainment.
