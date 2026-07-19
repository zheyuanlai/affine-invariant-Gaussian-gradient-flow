# Hilbert--Schmidt posterior curvature versus exact one-form coercivity

## 0. Verdict

Let \(Y=X+G\), where \(G\sim N(0,I_n)\), let

\[
 d\nu=e^{-U}dy,\qquad H=D^2U=I-\operatorname {Cov}(X\mid Y),
\]

and work in the intrinsic full-dimensional support.  Then

\[
 0\preceq H\preceq I,qquad
 \sum_j\left|(D_jH)u\right|^2\le16u^THu.           \tag{0.1}
\]

This note tests whether (0.1), by itself, proves the exact-gradient
coercivity

\[
 E|W|^2\le C\{E\|DW\|_{HS}^2+E[W^THW]\},
 \qquad W=\nabla f.                                \tag{OF}
\]

The answer from all direct multiplier and differentiated-eigenfunction
identities is negative: (0.1) removes the trace loss and gives
dimension-free local eigenspace coherence, but it does not supply the
global amplitude/nodal connectivity needed in (OF).  No counterexample is
claimed in the covariance-normalized Gaussian-output class; such a
counterexample would contradict KLS after reverse smoothing.  There is,
however, an exact Gaussian countermodel showing that the local differential
condition alone, without covariance normalization, cannot imply any
absolute coercivity constant.

## 1. Normalization and the exact energy

For the unscaled output, \(\operatorname {Cov}(Y)=2I\).  If
\(b=\nabla U\), integration by parts and matrix Cauchy--Schwarz give

\[
 EH=E[bb^T]\succeq\operatorname {Cov}(Y)^{-1}
 =\frac12I.                                        \tag{1.1}
\]

Thus a constant field is controlled.  For a gradient field, the integrated
Bochner identity is

\[
 \|(-L)f\|_2^2
 =\mathcal D(W)+\mathcal K(W),                     \tag{1.2}
\]

where

\[
 \mathcal D(W)=E\|DW\|_{HS}^2,qquad
 \mathcal K(W)=E[W^THW].                           \tag{1.3}
\]

Consequently (OF) is exactly the scalar Poincare inequality on the range
of the gradient.  The question is whether the additional local identity
(0.1) proves it.

## 2. What multiplication by \(H\) and \(I-H\) proves

Put \(P=I-H\).  From (0.1), \(\|H\|_{op},\|P\|_{op}\le1\), and
\(|a+b|^2\le2|a|^2+2|b|^2\), one gets for every \(H^1\) vector field

\[
 \begin{aligned}
 E\|D(HW)\|_{HS}^2&\le2\mathcal D+32\mathcal K,\\
 E\|D(PW)\|_{HS}^2&\le2\mathcal D+32\mathcal K,\\
 E|HW|^2&\le\mathcal K,\\
 \|W-PW\|_2^2&\le\mathcal K.                     \tag{2.1}
 \end{aligned}
\]

Hence a hypothetical normalized, almost-centered field with
\(\mathcal D+\mathcal K=o(1)\) produces another almost-normalized,
almost-centered field \(PW\) whose derivative is still \(o(1)\).  This is
a useful stability statement, but it reproduces the slow field rather than
controlling it.

The obstruction can be seen exactly for an eigen-gradient.  Suppose

\[
 (-L+H)W=\lambda W,qquad E|W|^2=1,                 \tag{2.2}
\]

so \(\mathcal D+\mathcal K=\lambda\).  Define

\[
 \begin{aligned}
 R_H&=E\sum_j(D_jW)^TH(D_jW),& R_P&=\mathcal D-R_H,\\
 S_H&=E[W^TH^2W],& S_P&=\mathcal K-S_H,\\
 \mathcal T&=E\sum_j\langle D_jW,(D_jH)W\rangle.
 \end{aligned}                                     \tag{2.3}
\]

Testing (2.2) weakly against \(HW\) and \(PW\), respectively, gives

\[
 \boxed{
 \lambda\mathcal K=R_H+S_H+\mathcal T,\qquad
 \lambda(1-\mathcal K)=R_P+S_P-\mathcal T.}        \tag{2.4}
\]

Their sum is only \(\lambda=\mathcal D+\mathcal K\).  The entire new
term is the cubic contraction \(\mathcal T\), and (0.1) gives exactly

\[
 \boxed{|\mathcal T|\le4\sqrt{\mathcal D\mathcal K}.}     \tag{2.5}
\]

Since \(4\sqrt{\mathcal D\mathcal K}\le
2(\mathcal D+\mathcal K)=2\lambda\), its constant is too large for
absorption in (2.4).  More generally, testing against
\((aI+bH)W\) only takes a linear combination of the two identities in
(2.4), so no affine matrix multiplier creates a new positive term.

## 3. What differentiating the eigen-equation adds

Let \(Z=DW=D^2f\).  For a smooth eigenfunction, differentiation of (2.2)
and full symmetry of \(D^3U\) give the exact tensor equation

\[
 -LZ+HZ+ZH+(D H)[W]=\lambda Z.                     \tag{3.1}
\]

Here \(((DH)[W])_{ij}=\sum_kU_{ijk}W_k\), and its inner product with
\(Z\) is exactly \(\mathcal T\).  Taking the \(L^2\) inner product with
\(Z\) yields

\[
 \boxed{
 \lambda\mathcal D
 =E\|DZ\|_{HS}^2+2R_H+\mathcal T.}                 \tag{3.2}
\]

Thus

\[
 \mathcal T\le\lambda\mathcal D,                 \tag{3.3}
\]

and eliminating \(\mathcal T\) between (2.4) and (3.2) gives

\[
 \boxed{
 E\|DZ\|_{HS}^2+R_H
 =\lambda(\mathcal D-\mathcal K)+S_H.}             \tag{3.4}
\]

Equation (3.3) controls the positive part of the multiplier cancellation
very strongly when \(\lambda\) is small.  It does not control the negative
part: (2.5) still permits \(\mathcal T=-c\lambda\) when
\(\mathcal D\) and \(\mathcal K\) are comparable fractions of
\(\lambda\).  Then (3.2) merely transfers order-\(\lambda\) energy to
the next derivative tensor.  All nonnegative scalar ledgers
(2.4)--(3.4) remain algebraically compatible with \(\lambda\downarrow0\).
Thus the differentiated equation does not yield a numerical gap without a
new sign, compactness, or global-coherence input.

The identities extend to exact eigenfunctions by elliptic regularity and
core approximation.  A bottom spectral-window vector has corresponding
form-error terms tending to zero with the window width; those errors do not
alter the obstruction above.

## 4. Local simple-eigenspace coherence and its exact limit

There is one substantial consequence of (0.1).  Suppose the lowest
eigenvalue \(h\) of \(H\) is simple and

\[
 \lambda_2(H)-h\ge\kappa>0,
\]

and let \(P_h\) be the rank-one lowest spectral projection.  Differentiating
the eigenprojection and using (0.1) gives

\[
 \sum_j\|D_jP_h\|_{HS}^2\le\frac{32h}{\kappa^2}.   \tag{4.1}
\]

For any \(W\), if \(Z_h=P_hW\), then

\[
 E|W-Z_h|^2\le\frac{\mathcal K}{\kappa},qquad
 E\|DZ_h\|_{HS}^2
 \le2\mathcal D+\frac{32}{\kappa^2}\mathcal K.   \tag{4.2}
\]

These estimates contain no factor of the dimension.  In the exact case
\(HW=0\), \(\ker H=\operatorname {span}\{W\}\), and \(W=\nabla f\),
(0.1) gives \((D_jH)W=0\).  Differentiating \(HW=0\) shows that every
column of \(DW\) lies in \(\operatorname {span}\{W\}\); symmetry of
\(DW\) forces

\[
 DW=c\,WW^T.                                       \tag{4.3}
\]

Hence the line of \(W\) is fixed on each connected nonvanishing component,
and its amplitude depends on one Euclidean coordinate.  This exact
rigidity is valid, but its quantitative version requires control across
regions where \(|W|\) is small.

More explicitly, under the probability

\[
 d\sigma=|W|^2d\nu/E|W|^2,
\]

the rank-one direction \(N=WW^T/|W|^2\) and the posterior projection
\(P_h\) have small mutual distance and small Dirichlet energies.  The
missing statement is a dimension-free global-coherence inequality forcing
them close to one fixed line.  The weight \(|W|^2\) is not log-concave and
may suppress the transition set, so neither log-concavity of \(\nu\) nor
(0.1) supplies that inequality.

Once a fixed line is obtained, curl-freeness makes the amplitude
one-dimensional, and the ordinary one-dimensional log-concave Poincare
inequality closes the branch.  Thus the precise missing hypothesis is
amplitude-weighted line-field connectivity (and its higher-rank analogue),
not a further local bound on \(DH\).

## 5. What the local condition alone cannot imply

Even actual Gaussian outputs show that (0.1) without covariance
normalization is insufficient.  Let

\[
 \nu_R=N(0,R^2I_n),\qquad R\ge1.
\]

It is the law of \(X+G\) with
\(X\sim N(0,(R^2-1)I_n)\).  Its Hessian is
\(H=R^{-2}I\), so \(DH=0\) and (0.1) holds with zero left side.  For the
constant gradient \(W=e_1=\nabla x_1\),

\[
 E|W|^2=1,qquad
 E\|DW\|_{HS}^2+E[W^THW]=R^{-2}.                  \tag{5.1}
\]

Thus no absolute one-form constant follows from the local differential
inequality and Gaussian-output representation alone.  Covariance
normalization, equivalently the average-curvature constraint (1.1), is
essential.

With \(\operatorname {Cov}(Y)=2I\) restored, no countermodel is produced
here.  In that normalized class, (OF) is the KLS-equivalent remaining
global assertion.  The audited Hilbert--Schmidt theorem has removed the
previous local trace obstruction, but a proof still needs a genuinely
global mechanism controlling amplitude-weighted changes of low-curvature
eigenspaces.
