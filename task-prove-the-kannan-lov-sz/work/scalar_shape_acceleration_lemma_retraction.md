# Retraction: the proposed scalar shape-acceleration inequality is false

The displayed inequality (1.2) in the scalar shape-acceleration note,
\[
m''(s)\,U_z(s,z)\le
U_{ss}+2m'U_{sz}+(m')^2U_{zz},
\]
does not follow from midpoint convexity and must not be used.

To see the algebra, write \(y=z-m(s)\) and compare
\(U(s+h,z)\), \(U(s-h,z)\) at fixed \(z\). The term involving
\(U_zm''h^2/2\) appears with the same sign in the midpoint expansion and
in the endpoint average, so it cancels. The resulting second-order
condition is only
\[
U_{ss}+2m'U_{sz}+(m')^2U_{zz}\ge0,
\]
which is the ordinary convexity condition along the affine path
\(s\mapsto(s,z+m'(s)s)\); it gives no separate bound on \(m''U_z\).

The pure-translation test makes the error explicit. For
\[
U(s,z)=\Phi(s)+W(z-m(s)),
\]
the affine-path curvature is
\[
A_s=\Phi''(s)-W'(z-m(s))\,m''(s).
\]
Joint convexity does not imply \(W'm''\le A_s\); the claimed inequality
would impose an unjustified extra factor of two. All conclusions in the
original scalar note that depend on (1.2)--(1.5), including its proposed
score-range bound for \(m''\), are retracted.

This retraction does not affect the independent smooth Gaussian-fiber
calculation in gaussian_full_wfi.md, whose Schur-complement argument uses
the exact matrix quadratic condition.
