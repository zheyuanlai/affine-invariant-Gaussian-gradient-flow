# Scalar score-range acceleration lemma

Let \(U\in C^2(J\times I)\) be jointly convex, with conditional density
\(q_s(z)=e^{-U(s,z)}/\rho(s)\) on an interval \(I_s\). Assume the conditional
mean \(m(s)\), variance \(\sigma_s^2\), and derivatives are finite. Put
\(b=m'(s)\) and
\[
 A_s(z)=U_{ss}(s,z)+2b\,U_{sz}(s,z)+b^2U_{zz}(s,z)\ge0.       \tag{1.1}
\]

## Midpoint inequality

For \(h>0\), use the two points
\((s-h,m(s-h)+y)\) and \((s+h,m(s+h)+y)\). Their midpoint is
\((s,m(s)+y+\delta_h)\), with
\(\delta_h=\frac12[m(s-h)+m(s+h)-2m(s)]\). Joint convexity and Taylor
expansion give
\[
 U_z(s,m(s)+y)\,m''(s)\le A_s(m(s)+y).                       \tag{1.2}
\]
(The inequality is distributional for BV \(m'\), by finite differences.)

Let \(z_0\) be a mode of \(q_s\). If \(m''(s)\ge0\), integrate (1.2) over
the right side of the mode. Since \(U_z\) is nondecreasing and
\(q_s\propto e^{-U}\),
\[
 E_s[U_z\,{\bf1}_{\{z\ge z_0\}}]
 =q_s(z_0)-\lim_{z\to\sup I_s}q_s(z)\ge q_s(z_0).          \tag{1.3}
\]
For \(m''<0\), integrate \(-U_z\) on the left side and use the analogous
identity. A centered one-dimensional logconcave density of variance
\(\sigma_s^2\) has \(q_s(z_0)\ge3/(16\sigma_s)\): Chebyshev puts mass at
least \(3/4\) in an interval of length \(4\sigma_s\), so its maximum is at
least \(3/(16\sigma_s)\). Consequently
\[
 \boxed{\ |m''(s)|\le \frac{16}{3}\,\sigma_s\,E_s A_s.\ }   \tag{1.4}
\]

Using the affine-shear identity (shape_centroid_affine_shear_addendum.md),
\[
 E_sA_s=\Phi''(s)+C_s+B_s^0,                              \tag{1.5}
\]
where \(F_s^0=F_s+m'(s)\), \(C_s=E_s|D_zF_s^0|^2\), and
\(B_s^0=E_s[U_{zz}(F_s^0)^2]\). Thus scalar shape acceleration is controlled
by the marginal curvature, deformation charge, and exactly one residual
curvature term. No conditional Poincare inequality is used.

The factor \(\sigma_s\) and \(B_s^0\) are genuine remaining obstructions:
the cone stress test in work/gaussian_curvature_korn.md has
\(\sigma_s\asymp a\) and \(B_s^0\asymp a^{-1}\) in its boundary layer.
Therefore (1.4) alone does not yield the desired integrated WFI bound.

