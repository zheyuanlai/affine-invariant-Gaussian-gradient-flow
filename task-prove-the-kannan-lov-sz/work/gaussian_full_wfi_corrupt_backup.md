# Full weighted-Fisher bound for smooth Gaussian conditional fibers

This note records a complete calculation for the Gaussian conditional-slice
subclass.  It is dimension-free and allows (R'(s)) not to commute with
(R(s)).  The lower covariance bound in the statement is the one supplied
by unit transverse Gaussian noise; without it, the centroid term below is
not controlled by the hypotheses used here.

## Statement

Let

\[
 p(s,z)=Z^{-1}\exp\left\{-W(s)-\frac12(z-m(s))^TQ(s)(z-m(s))\right\},
 \qquad R(s)=Q(s)^{-1}.
\]

Assume (W,m,R) are (C^2), (R(s)succ0), the joint potential is convex,
and the one-dimensional marginal (ho=e^{-\Phi}) is centered with variance
one.  Let (\tau) be its canonical Stein kernel.  Assume

\[
 I\preceq R(s)\quad\text{for all }s,
 \qquad E R(S)\preceq\Lambda I,
 \qquad E[S,m(S)]=0.
\]

For the centered conditional Poisson field (L_sg_s=\partial_s\log q_s),
write (F_s=\nabla_zg_s),
\[
 C_s=E_s\|D_zF_s\|_{\mathrm{HS}}^2,
 \qquad B_s=E_s\langle Q(s)F_s,F_s\rangle .
\]
Then

\[
 \int \rho(s)\tau(s)^2\bigl(C_s+B_s\bigr)\,ds
 \le \frac32+C_\mathrm{cent}\Lambda,
\]

where (C_\mathrm{cent}) is the universal constant from the Gaussian
centroid estimate (one may take the explicit, deliberately loose value
appearing in `shape_centroid_stability_clean.md`).

The estimate is unchanged if (R,R') are only distributionally defined and
the displayed quantities are obtained by a smooth approximation with the same
pointwise bounds; the approximation issue for a general log-concave law is
separate.

## Exact Poisson field

Put (y=z-m(s)).  The conditional score is

\[
 \ell_s=m'(s)^TQy-rac12\left(y^TQ'y-\operatorname{tr}(Q'R)\right).
\]

If (A=A^T) is defined by

\[
 QA+AQ=Q'=-QR'Q,
 \tag{2.1}
\]

then (g_s=-m'(s)^Tz+\frac12y^TAy) (up to an additive constant) solves
\(L_sg_s=\ell_s\).  Consequently

\[
 F_s=-m'(s)+Ay,
 \qquad C_s=\|A\|_{\mathrm{HS}}^2,
 \tag{2.2}
\]

and, since (E_sy=0),

\[
 B_s=m'(s)^TQm'(s)+B_s^0,
 \qquad B_s^0:=\operatorname{tr}(QARA).
 \tag{2.3}
\]

## Noncommuting matrix estimate

Diagonalize (R) at the fixed (s): (R=\operatorname{diag}(r_1,\ldots,r_d)).
No derivative of the diagonalizing basis is taken.  Equation (2.1) becomes

\[
 A_{ij}=-\frac{R'_{ij}}{r_i+r_j}.
 \tag{2.4}
\]

Define the covariance Fisher term

\[
 J:=\frac12\operatorname{tr}(R^{-1}R'R^{-1}R')
   =\frac12\sum_{i,j}\frac{(R'_{ij})^2}{r_ir_j}.
 \tag{2.5}
\]

From (2.3)--(2.4),

\[
 B_s^0=\sum_{i,j}\frac{r_j}{r_i}\frac{(R'_{ij})^2}{(r_i+r_j)^2}.
 \tag{2.6}
\]

For (i=j), the summand in (2.6) is one half of the corresponding
diagonal contribution to (J).  For each unordered pair (i<j), combining
the two ordered summands gives

\[
 (R'_{ij})^2\frac{r_i^2+r_j^2}{r_ir_j(r_i+r_j)^2}
 \le \frac{(R'_{ij})^2}{r_ir_j},
\]

which is exactly the pair contribution to (J).  Therefore

\[
 0\le B_s^0\le J.                                      \tag{2.7}
\]

Likewise (pairing ordered terms in (C_s=\sum A_{ij}^2)) gives

\[
 0\le C_s\le \frac12J.                                 \tag{2.8}
\]

Joint convexity gives the Schur-complement formula

\[
 \Phi''=W''+\frac12\operatorname{tr}(R^{-1}(-R''))+J\ge J.
 \tag{2.9}
\]

Thus pointwise

\[
 C_s+B_s^0\le\frac32\Phi''(s).                         \tag{2.10}
\]

The canonical Stein identity for centered variance-one log-concave (\rho)
is

\[
 \int\rho\tau^2\Phi''=1-E_\rho(\tau')^2\le1.          \tag{2.11}
\]

Hence

\[
 \int\rho\tau^2(C_s+B_s^0)\,ds\le\frac32.             \tag{2.12}
\]

## Centroid term

The mixed-isotropy identity (E[\tau(S)m'(S)]=0), the Schur complement, and
the concavity of (R) yield the dimension-free estimate proved in
`shape_centroid_stability_audited.md`:

\[
 E[\tau(S)^2|m'(S)|^2]\le C_\mathrm{cent}\Lambda .       \tag{2.13}
\]

The unit-noise lower bound (R\succeq I) implies (Q\preceq I), so

\[
 \int\rho\tau^2,m'^TQm'\,ds
 \le E[\tau^2|m'|^2]\le C_\mathrm{cent}\Lambda .          \tag{2.14}
\]

Combining (2.3), (2.12), and (2.14) proves the statement.

## Scope warning

The argument is a subclass theorem.  For an arbitrary Gaussian family with
only (ER\preceq\Lambda I), but no uniform lower bound on (R), (2.14) is
unavailable: (m'^TQm') can be much larger than (|m'|^2).  Unit Gaussian
convolution is therefore part of the hypotheses when this lemma is used in a
post-noise reduction.  Extending the estimate uniformly through removal of
that noise, or to non-Gaussian conditional fibers, remains the unresolved
KLS-strength step.
