# Checkpoint 15 patch text: Gaussian-fiber promotion and scalar-note retraction

The following paragraphs should be inserted into checkpoint_15.md.

## Candidate-proof status replacement

There is still no complete proof of KLS, and no KLS conclusion is claimed at
this checkpoint. The new work proves a dimension-free forward
self-convolution theorem, closes the weighted-Fisher estimate for every pure
translating slice family, and proves a full weighted-Fisher estimate for
smooth full-support Gaussian conditional fibers after unit transverse
noise. It identifies two exact global gates that survive the present audits
and supplies explicit counterexamples to three tempting reverse or
slice-local arguments.

## Addition after the general weighted ledger

For the smooth Gaussian conditional subclass \(q_s=N(m(s),R(s))\), with
\(R(s)\succeq I\), the residual can be evaluated exactly even when \(R'\) and
\(R\) do not commute. Writing \(F_s=-m'(s)+A_sy\), with
\(A_sR+RA_s=-R'\), and
\[
J_s=\frac12\operatorname{tr}(R^{-1}R'R^{-1}R'),
\]
diagonalization of \(R(s)\) and pairing the \(i,j\) terms gives
\[
C_s+B_s^0=J_s,\qquad
B_s=m'(s)^TR(s)^{-1}m'(s)+B_s^0.
\]
The marginal curvature is
\[
\Phi''=W''+\frac12\operatorname{tr}(R^{-1}(-R''))+J_s\ge J_s.
\]
Therefore \(\int\rho\tau^2(C_s+B_s^0)\,ds\le1\). The audited centroid
estimate \(E[\tau^2|m'|^2]\le C\Lambda\), together with \(R\succeq I\),
controls the remaining term \(m'^TR^{-1}m'\), and
\[
\int\rho\tau^2(C_s+B_s)\,ds\le1+C\Lambda.
\]
This is a \(C^2\), full-support, unit-noise subclass theorem only; no
nonsmooth or zero-noise limit is asserted.

## Shape-centroid section replacement

Replace the sentence “No nonsmooth limit or full \(B_s\) estimate is claimed
here” in item 1 by:

In the post-noise normalization \(R\succeq I\), the exact Gaussian Poisson
calculation gives
\[
C_s+B_s^0
=\frac12\operatorname{tr}(R^{-1}R'R^{-1}R')
\le\Phi''(s),\qquad
B_s=m'^TR^{-1}m'+B_s^0.
\]
Consequently
\[
\int\rho\tau^2(C_s+B_s)\,ds\le1+C\Lambda.
\]
This allows noncommuting, rotating covariance matrices. It is a \(C^2\),
full-support, unit-noise subclass result; no nonsmooth or zero-noise limit
is claimed.

## Mandatory retraction

The proposed scalar shape-acceleration inequality in
scalar_shape_acceleration_lemma.md is false. Midpoint expansion cancels the
\(U_zm''\) term and yields only affine-path curvature \(A_s\ge0\). See
scalar_shape_acceleration_lemma_retraction.md; all uses of equations
(1.2)--(1.5) in the original scalar note must be removed.
