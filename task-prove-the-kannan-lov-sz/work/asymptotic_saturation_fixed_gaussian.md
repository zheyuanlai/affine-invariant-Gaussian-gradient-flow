# Asymptotic saturation of the fixed-Gaussian comparison

## 0. Statement

Let \(\mu_n\) be centered isotropic log-concave laws, let
\(a_n=\lambda_1(\mu_n)\), let \(X_n\sim\mu_n\), \(G_n\sim N(0,I)\), and
write
\[
 Y_n=X_n+G_n,\qquad \lambda_n=\lambda_1(\mathcal L(Y_n)).
\]
The normalized output \(S_n=Y_n/\sqrt2\) has gap \(b_n=2\lambda_n\).
The posterior reverse comparison and the first-order posterior lift imply,
for every \(a_n<1\),
\[
 \frac{\sqrt{1+12a_n}-1}{6}\ \le\ \lambda_n\
 \le \frac{a_n}{1-a_n}.                              \tag{0.1}
\]
Consequently, if \(a_n\to0\),
\[
 a_n-3a_n^2\le\lambda_n\le a_n+2a_n^2
 \quad\text{for all sufficiently large \(n\)},       \tag{0.2}
\]
and hence
\[
 \lambda_n=a_n+O(a_n^2),\qquad b_n=2a_n+O(a_n^2).      \tag{0.3}
\]
The constants in (0.1)--(0.3) are dimension free.

The same comparison forces near-saturation of the posterior conditional
Poincare and covariance inequalities for an input bottom-window vector.
It does **not**, by itself, force a common posterior Gaussian factor or a
positive ANOVA interaction: the input and output minimizing vectors need
not coincide, and the deficits are only averaged, with no control on rare
posterior fibers.

## 1. The two inequalities

### 1.1 Reverse posterior comparison

For a centered unit input test \(g\) with energy \(q<1\), put
\(F(y)=E[g(X)\mid Y=y]\), and
\[
 d=E\operatorname {Var}(g(X)\mid Y),\qquad
 e=E|\nabla F(Y)|^2.
\]
Posterior strong log-concavity and covariance Cauchy--Schwarz give
\[
 d\le q,\qquad e\le d,\qquad
 \operatorname {Var}(F)=1-d.                         \tag{1.1}
\]
The output Rayleigh quotient is at most \(e/(1-d)\), so, at the spectral
edge,
\[
 \lambda_n\le {q\over1-q}.
\]
Taking \(q\downarrow a_n\) and rearranging gives
\[
 a_n\ge{\lambda_n\over1+\lambda_n},
 \qquad\lambda_n\le {a_n\over1-a_n}.                 \tag{1.2}
\]
This passage is valid for non-attained edges by spectral-window
approximation.

### 1.2 First-order posterior lift

Let \(f\) be a centered unit vector in a bottom spectral window of the
unscaled output, and set
\[
 q=E|\nabla f|^2,\quad W=\nabla f,\quad
 B=E\|D^2f\|_{\rm HS}^2,\quad C=E[W^T(I-A_Y)W].
\]
For an exact eigenvector with eigenvalue \(\lambda_n\),
\[
 q=\lambda_n,\qquad B+C=\lambda_n^2.                 \tag{1.3}
\]
The conditional lift
\[
 \Phi=f(Y)+W(Y)\cdot(X-m(Y)),\qquad
 u(X)=E_G[\Phi\mid X]
\]
satisfies
\[
 \|u\|_2^2\ge1+q-2(B+C),\qquad
 E|\nabla u|^2\le q+3(B+C).                           \tag{1.4}
\]
Therefore
\[
 a_n\le Q(\lambda_n):=
 \frac{\lambda_n+3\lambda_n^2}
      {1+\lambda_n-2\lambda_n^2}.                    \tag{1.5}
\]
For a spectral window, replace \(\lambda_n\) on the right by
\(b+\varepsilon\) and pass \(\varepsilon\downarrow0\).

For \(0\le\lambda\le1/2\),
\[
 Q(\lambda)\le\lambda+3\lambda^2.                    \tag{1.6}
\]
Combining (1.2) and (1.6) gives \(a_n\le\lambda_n+3\lambda_n^2\), hence
\[
 \lambda_n\ge {\sqrt{1+12a_n}-1\over6}.               \tag{1.7}
\]
The elementary bounds
\[
 a-3a^2\le{\sqrt{1+12a}-1\over6},\qquad
 {a\over1-a}\le a+2a^2\quad(0\le a\le1/2)
\]
give (0.2).

## 2. Deficit ledger for the reverse transfer

Choose a centered unit \(g_n\) with
\(q_n\le a_n+o(a_n^2)\); this is possible by taking a spectral-window
vector with window width \(o(a_n^2)\). Let \(d_n,e_n\) be as in (1.1).
Set
\[
 R_n={q_n\over1-q_n}-\lambda_n.
\]
By (0.2),
\[
 R_n\le 4a_n^2+o(a_n^2).                              \tag{2.1}
\]
The monotone chain
\[
 \lambda_n\le {e_n\over1-d_n}\le {d_n\over1-d_n}
 \le {q_n\over1-q_n}
\]
then gives
\[
 q_n-d_n=O(a_n^2),\qquad d_n-e_n=O(a_n^2),\qquad
 {e_n\over1-d_n}-\lambda_n=O(a_n^2).                  \tag{2.2}
\]
More explicitly, for \(a_n\le1/4\) each quantity is bounded by
\(8a_n^2+o(a_n^2)\).

Write \(\pi_y=\mathcal L(X\mid Y=y)\), \(A_y=\operatorname{Cov}_{\pi_y}X\),
and \(c_y=\operatorname{Cov}_{\pi_y}(X,g)\).  Then
\[
 e_n=E|c_Y|^2,\qquad d_n=E\operatorname{Var}_{\pi_Y}g.
\]
The two nonnegative conditional deficits split as
\[
\begin{aligned}
q_n-d_n
 &=E\bigl[E_{\pi_Y}|\nabla g|^2-\operatorname{Var}_{\pi_Y}g\bigr],\\
d_n-e_n
 &=E\bigl[\operatorname{Var}_{\pi_Y}g-|c_Y|^2\bigr].  \tag{2.3}
\end{aligned}
\]
Thus both the conditional \(1\)-strongly-log-concave Poincare deficit and
the covariance Cauchy--Schwarz deficit are \(O(a_n^2)\), while
\(d_n,e_n=a_n+O(a_n^2)\).

A more geometric nonnegative split of the second line is
\[
\operatorname{Var}_{\pi_y}g-|c_y|^2
=
\bigl(\operatorname{Var}_{\pi_y}g
      -c_y^TA_y^\dagger c_y\bigr)
+
c_y^T(A_y^\dagger-I)c_y,                              \tag{2.4}
\]
on the support of \(A_y\). The first term is the posterior covariance
Cauchy--Schwarz/BL deficit and the second measures failure of \(c_y\) to
lie in eigen-directions of \(A_y\) with eigenvalue \(1\). Their integrals
are each \(O(a_n^2)\). In particular, the posterior signal carried by
\(g_n\) is concentrated (in the weighted \(L^2\) sense) on nearly
saturating covariance directions.

## 3. Deficit ledger for the output lift

For an output bottom-window vector \(f_n\), let
\[
 \rho_n=\|{\cal A}f_n\|_2^2,\qquad
 B_n=E\|D^2f_n\|_{\rm HS}^2,\qquad
 C_n=E[\nabla f_n^T(I-A_Y)\nabla f_n].
\]
Then \(\rho_n=B_n+C_n=\lambda_n^2+o(\lambda_n^2)\). The lift
\(\Phi_n=f_n+\nabla f_n\cdot(X-m)\) has the exact identities
\[
 E[\Phi_n\mid Y]=f_n,\qquad
 \|\Phi_n\|_2^2=1+q_n-C_n,                           \tag{3.1}
\]
and
\[
 E|\nabla_G\Phi_n|^2
 =E|H\nabla f_n+D^2f_n\,Z|^2
 \le \rho_n.                                        \tag{3.2}
\]
The loss in conditional Gaussian Poincare,
\[
 \Delta_{G,n}:=\rho_n-
 E|\nabla_G\Phi_n|^2\ge0,                              \tag{3.3}
\]
is not controlled by (0.1); it includes the slack in \(H^2\preceq H\) and
\(A\preceq I\). The input projection \(u_n=E_G\Phi_n\) obeys
\[
 \|u_n\|_2^2
 =1+q_n-C_n-\operatorname{Var}(\Phi_n\mid X),\qquad
 E|\nabla u_n|^2\le q_n+3\rho_n.                      \tag{3.4}
\]
Since \(\lambda_n=a_n+O(a_n^2)\), its Rayleigh quotient is
\[
 {E|\nabla u_n|^2\over\|u_n\|_2^2}
 \le a_n+O(a_n^2).                                   \tag{3.5}
\]
Consequently the lift gives another near-bottom input vector, but (3.5)
does not identify it with the reverse-transfer vector \(g_n\). No estimate
in the two chains forces
\[
 \operatorname{Var}(\Phi_n\mid X)=O(a_n^2)
 \quad\text{or}\quad
 C_n=o(a_n^2);
\]
only the upper bounds \(O(a_n^2)\) are available.

## 4. Why simultaneous saturation does not yet force a Gaussian factor

Equality in the conditional \(1\)-Poincare inequality for a posterior
\(\pi_y\) would imply a Gaussian factor in the selected direction (and an
affine conditional profile). However, (2.3)--(2.4) provide only an
integrated \(O(a_n^2)\) deficit while the total conditional variance is
\(O(a_n)\). They allow an \(O(a_n)\) fraction of posterior mass to carry
order-one relative deficit, and they do not synchronize the selected
directions \(c_y\) across \(y\).

The HS saturation theorem controls this synchronization only in the
following weighted sense:
\[
 \sum_j\|(DA_y[e_j])u\|^2
 \le16\,u^T(I-A_y)u.                                  \tag{4.1}
\]
For the selected \(u=c_y\)-direction, integration of (4.1) is
\(O(a_n^2)\) after the natural normalization. It does not bound the
variation of the direction on fibers where the posterior signal is small,
nor does it rule out rare switching between different nearly-saturating
directions. A quantitative equality/stability theorem converting
(2.3)--(4.1) into a global Gaussian factor would itself be a new
dimension-free rigidity result; none follows from the present estimates.

Likewise, the output lift produces a near-bottom input vector but not an
identity between its Gaussian ANOVA residual and that of \(g_n\). Therefore
the simultaneous \(O(a_n^2)\) scalar saturation does not presently imply
\[
 E R_n^2\ge c\,\operatorname{dist}(f_n,\mathrm{Aff})^2
\]
for the fixed-Gaussian ANOVA residual. Establishing such a lower bound is
exactly the unresolved interaction gate.

## 5. Formal edge cases

The argument uses intrinsic Gaussian noise on the supporting affine hull.
Spectral-window vectors avoid assuming an attained gap. If \(a_n=0\) for a
non-point probability on a fixed finite-dimensional support, the usual
one-dimensional/log-concave Cheeger positivity excludes this; the asymptotic
statement concerns \(a_n\downarrow0\) across dimensions. All limits above
are taken first in the spectral-window width and then along the sequence,
with constants independent of the dimension.

