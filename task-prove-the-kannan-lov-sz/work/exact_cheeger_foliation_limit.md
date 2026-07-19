# Exact Euclidean Cheeger foliation obtained from vanishing deficit

## 0. Purpose and verdict

The constant-anisotropy replacement is needed only at a positive error
scale. If the heat construction is run along a sequence for which its
coarea Cheeger deficit tends to zero, compactness removes the anisotropy
altogether. The abstract theorem below gives an exact Euclidean limiting
function whose complete coarea normal matrix is the prescribed limiting
matrix.

This is a strict strengthening of merely producing near-Cheeger leaves:
almost every nontrivial level of the limit is an attained Cheeger set. It
does not solve the remaining balanced one-interface branch. A two-valued
function supported on a balanced Cheeger set is an exact example of that
branch.

## 1. Setting

Let \(\mu\) be a fixed non-atomic log-concave probability on its
\(k\)-dimensional affine support, \(k\geq1\), and let
\(p=\psi_\mu>0\). Weighted variation and perimeter are relative to the
affine support. For \(U\in BV(\mu)\), set

\[
 R(U)=\inf_{c\in\mathbb R}\int |U-c|\,d\mu,
 \qquad D_p(U)={\rm TV}_\mu(U)-pR(U)\geq0,                 \tag{1.1}
\]

and, with \(DU=\sigma_U|DU|_\mu\),

\[
 M(U)=\int\sigma_U\sigma_U^T\,d|DU|_\mu,
 \qquad {\rm tr}\,M(U)={\rm TV}_\mu(U).                  \tag{1.2}
\]

All functions are considered modulo additive constants and represented by
a median-zero version. Then \(R(U)=\|U\|_{L^1(\mu)}\).

## 2. The compactness theorem

**Theorem 2.1 (vanishing-deficit matrix retention).** Suppose that
\(F_j\in BV(\mu)\) and a positive semidefinite matrix \(M_0\) satisfy

\[
 D_p(F_j)\longrightarrow0,
 \qquad M(F_j)\longrightarrow M_0
       \quad\hbox{in nuclear norm},                         \tag{2.1}
\]

where \(0<{\rm tr}\,M_0<\infty\). Fix any
\(0<\kappa<1/3\). For every \(j\), let \(G_j\) satisfy

\[
 D_p(G_j)+\kappa\|M(G_j)-M(F_j)\|_*
       \le D_p(F_j)+j^{-1}D_p(F_j).                         \tag{2.2}
\]

Such \(G_j\) may be obtained either from the exact boxed minimization after
putting both functions in a sufficiently large common value interval, or
from the free-law Ekeland construction. Assume in addition that, after
choosing median-zero representatives, \((G_j)\) is uniformly integrable in
\(L^1(\mu)\). This holds in particular when all \(G_j\) lie in one fixed
value interval. Then, after passing to a
subsequence and adding constants, there is \(G\in BV(\mu)\) such that

\[
 G_j\to G\quad\hbox{in }L^1(\mu),
 \qquad D_p(G)=0,
 \qquad M(G)=M_0.                                         \tag{2.3}
\]

In particular the convergence is strict in weighted \(BV\):

\[
 {\rm TV}_\mu(G_j)\longrightarrow{\rm TV}_\mu(G)
       ={\rm tr}\,M_0.                                    \tag{2.4}
\]

**Proof.** Put \(\varepsilon_j=(1+j^{-1})D_p(F_j)\). The two
nonnegative terms in (2.2) give

\[
 D_p(G_j)\le\varepsilon_j\to0,
 \qquad
 \|M(G_j)-M(F_j)\|_*\le\varepsilon_j/\kappa\to0.         \tag{2.5}
\]

Taking traces in the second assertion and using (2.1),

\[
 {\rm TV}_\mu(G_j)={\rm tr}\,M(G_j)\longrightarrow
                         {\rm tr}\,M_0.                    \tag{2.6}
\]

Choose a median-zero representative. By the definition of \(p\),

\[
 \|G_j\|_1=R(G_j)\le p^{-1}{\rm TV}_\mu(G_j),             \tag{2.7}
\]

so both the quotient \(L^1\) norm and the variation are uniformly bounded.
On every compact subset of the relative interior, the log-concave density
is bounded above and below by positive constants. Local \(BV\)
compactness and a diagonal extraction give local convergence in measure.
The assumed uniform integrability, together with tightness of \(\mu\), then
upgrades this convergence to global \(L^1(\mu)\) convergence. Thus,
after a further choice of constants if necessary,
\(G_j\to G\) in \(L^1(\mu)\).

The quotient norm \(R\) is one-Lipschitz in \(L^1\), hence
\(R(G_j)\to R(G)\). From (1.1), (2.5), and (2.6),

\[
 pR(G)=\lim_j\bigl({\rm TV}_\mu(G_j)-D_p(G_j)\bigr)
      ={\rm tr}\,M_0.                                     \tag{2.8}
\]

Lower semicontinuity gives
\({\rm TV}_\mu(G)\le{\rm tr}\,M_0\), whereas the Cheeger inequality
gives \({\rm TV}_\mu(G)\ge pR(G)={\rm tr}\,M_0\).
Thus equality holds, proving (2.4) and \(D_p(G)=0\).

Strict weighted-\(BV\) convergence and Reshetnyak continuity apply to every
continuous function of the polar direction. Entrywise, with
\(h_{ab}(n)=n_an_b\), this yields

\[
 M(G_j)\longrightarrow M(G).                              \tag{2.9}
\]

Equations (2.1) and (2.5) identify the same limit as \(M_0\). This proves
(2.3). \(\square\)

## 3. Exact levelwise consequence

**Corollary 3.1 (exact Cheeger foliation).** Under Theorem 2.1, for almost
every \(r\in\mathbb R\), writing \(A_r=\{G>r\}\),

\[
 P_\mu(A_r)=p\min\{\mu(A_r),1-\mu(A_r)\}.                  \tag{3.1}
\]

Moreover,

\[
 M_0=M(G)=\int_{\mathbb R}M(A_r)\,dr.                     \tag{3.2}
\]

**Proof.** Weighted coarea and layer cake give

\[
 0=D_p(G)=\int_{\mathbb R}
 \left[P_\mu(A_r)-p\min\{\mu(A_r),1-\mu(A_r)\}\right]dr. \tag{3.3}
\]

The integrand is nonnegative for every finite-perimeter level, so it
vanishes almost everywhere. Matrix-valued coarea gives (3.2). \(\square\)

If \(Q_0=M_0/{\rm tr}M_0\), every lower bound on
\(1-{\rm tr}(Q_0^2)\) is therefore retained with no anisotropic conversion
loss and no singular-ray approximation.

## 4. Diffuse levels versus the balanced plateau

Assume temporarily that the density is smooth and positive on all of the
affine support. If a level in (3.1) has volume different from \(1/2\),
then its volume remains on the same linear branch of the Cheeger tent under
all sufficiently small compactly supported deformations. Consequently it
is an unconstrained local minimizer of

\[
 P_\mu(B)-p\mu(B)\quad(\mu(A_r)<1/2),
 \qquad
 P_\mu(B)+p\mu(B)\quad(\mu(A_r)>1/2).                      \tag{4.1}
\]

On a regular component, the weighted Jacobi inequality holds for every
compactly supported normal lapse. The cutoff argument from the direct
deficit audit then forces

\[
 S=0,
 \qquad \nabla^2V[n,n]=0                                  \tag{4.2}
\]

almost everywhere on that component. Thus every genuinely diffuse exact
level is flat and log-affine in its normal direction.

The exceptional set of values

\[
 I_{1/2}=\{r:\mu(G>r)=1/2\}                                \tag{4.3}
\]

is an interval. For any two of its interior values, their superlevel sets
agree modulo \(\mu\), since the measure of their set difference is the
difference of their distribution values. Hence all variation contributed
by the interior of (4.3) is a multiple of the normal matrix of one balanced
Cheeger interface. A two-valued function on such an interface shows that
this contribution may equal the whole matrix in (3.2).

Therefore the limiting theorem reduces the geometric problem exactly to:

1. classify the flat/log-affine diffuse components (rank-one in the
   smooth full-support case); and
2. prove a dimension-free inverse for one balanced Cheeger interface whose
   normalized normal matrix has universal angular variance.

The second item remains load bearing.

## 5. How the heat sequence is to be instantiated

The fixed-scale heat estimates have the schematic form

\[
 {D_p(F_\alpha)\over p}
   \le \beta_\alpha+C\sqrt\alpha,
 \qquad
 {{\rm tr}M(F_\alpha)\over p}
   \ge {c\over\sqrt{\log(e/\alpha)}},                      \tag{5.1}
\]

and the normalized matrices retain a universal angular-variance lower
bound in the high-\(C_P\) branch. Multiplying \(F_\alpha\) by

\[
 a_\alpha={p\over{\rm tr}M(F_\alpha)}                      \tag{5.2}
\]

normalizes its matrix trace to \(p\), while one-homogeneity gives

\[
 {D_p(a_\alpha F_\alpha)\over p}
 ={D_p(F_\alpha)\over{\rm tr}M(F_\alpha)}
 \le C\sqrt{\alpha\log(e/\alpha)}+o(1).                  \tag{5.3}
\]

Thus the deficit tends to zero without losing the normalized normal law.
To invoke Theorem 2.1 one still has to extract a nuclear-norm limit of those
normalized matrices and verify uniform integrability of the scaled
comparators/replacements. In a fixed finite dimension matrix compactness is
ordinary, but uniform integrability is a separate requirement: bounded
quotient \(L^1\) norm alone does not prevent value spikes escaping into the
spatial tails.
In a contradiction sequence with growing dimension, (5.1)--(5.3) are used
measure by measure; they do not by themselves identify a common ambient
matrix across dimensions.

The remaining high-rank inverse cannot be hidden in this limiting passage:
its constants must be uniform in the dimension of each individual measure.
