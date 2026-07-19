# Affine-shear identity for the centroid sector

This addendum fixes the signs in the affine-shear calculation used in the
shape-centroid audit. It is a local identity for a regular conditional
slice; it is not a KLS estimate.

Assume \(U\in C^2\), \(q_s(z)=e^{-U(s,z)}/\rho(s)\), \(\rho=e^{-\Phi}\), with
enough decay (or boundary conditions) for integration by parts. Put
\(H_s=D^2_{zz}U\), \(L_s=\Delta_z-\nabla_zU\cdot\nabla_z\), and let
\(g_s\) solve
\[
 L_sg_s=\ell_s,\qquad \ell_s=-U_s+\Phi',\qquad F_s=\nabla_zg_s.
\]
Then \(E_s\nabla_zU=0\), \(E_sF_s=-m'(s)\). Set \(b=m'(s)\) and
\(F_s^0=F_s+b\), so \(E_sF_s^0=0\).

## Identity

Use the fixed affine shear \(z=y+bs\), and define
\[
 \widetilde U_b(s,y)=U(s,y+bs).
\]
Its \(s\)-derivatives at the chosen value of \(s\) are
\[
 \partial_s\widetilde U_b=U_s+b\cdot U_z,\qquad
 \partial_{ss}\widetilde U_b
 =U_{ss}+2U_{sz}\cdot b+b^TH_sb=:A_b.                 \tag{1.1}
\]
The shear has unit Jacobian, so it leaves \(\rho\) unchanged. Differentiating
\(\rho(s)=\int e^{-\widetilde U_b(s,y)}dy\) gives
\[
 \Phi''=E_sA_b-\operatorname{Var}_s(U_s+b\cdot U_z).    \tag{1.2}
\]
Since \(E_sU_z=0\) and \(E_sU_s=\Phi'\),
\[
 \operatorname{Var}_s(U_s+b\cdot U_z)
 =E_s(\ell_s-b\cdot U_z)^2.                              \tag{1.3}
\]
Moreover \(L_s(b\cdot z)=-b\cdot U_z\), hence
\[
 L_s(g_s+b\cdot z)=\ell_s-b\cdot U_z,\qquad
 \nabla_z(g_s+b\cdot z)=F_s^0.                           \tag{1.4}
\]
The integrated Bochner identity therefore yields
\[
 E_s(\ell_s-b\cdot U_z)^2
 =E_s\|D_zF_s^0\|_{\rm HS}^2
  +E_s\langle H_sF_s^0,F_s^0\rangle
 =C_s+B_s^0.                                             \tag{1.5}
\]
Combining (1.2)--(1.5) gives the exact affine-shear formula
\[
 \boxed{\quad
 E_s\,D^2U[(1,m'(s)),(1,m'(s))]
 =\Phi''(s)+C_s+B_s^0,\quad}                             \tag{1.6}
\]
where \(B_s^0=E_s\langle H_sF_s^0,F_s^0\rangle\).

## Consequence and limitation

Formula (1.6) is the precise algebraic split of the centroid and shape
sectors. It rules out replacing the affine-shear curvature mass by
\(\Phi''+C_s\) without separately controlling \(B_s^0\). It does **not**
show that a global weighted inequality is impossible: \(B_s^0\) can be
charged by joint \(s\)-geometry and the Stein weight. The cone example in
work/gaussian_curvature_korn.md has \(m'=0\), \(B_s^0=B_s\) large on a
boundary layer, yet its full weighted integral is bounded. Thus (1.6)
only identifies the exact remaining term; any claim that the deformation
budget alone controls it would be an unproved additional lemma.

