# Checkpoint 14: eigenfunction graph rigidity and slice-curvature compensation

## 1. Candidate-proof status

There is not yet a complete dimension-free proof, and no KLS conclusion is
claimed at this checkpoint.  The new work isolates one narrow second-order
gate for a near-linear first eigenfunction, proves that gate on several
genuinely dependent classes, and shows exactly why ordinary thin-shell
information stops at a logarithm.

Let \(\mu\) be isotropic and log-concave, and let

\[
 -L_\mu f=\lambda f,\qquad
 \mathbb Ef=0,\qquad \mathbb Ef^2=1,\qquad
 \mathbb E|\nabla f|^2=\lambda .
\]

Write

\[
 a=\mathbb E[Xf(X)],\qquad
 \delta^2=1-|a|^2,\qquad r=f-a\cdot x .
\]

The tensor-extremizer route would make \(\delta\) numerical and small if its
convex-potential normal-cone term could be eliminated.  The present
checkpoint studies the independent question of whether such near-linearity
already forces \(\lambda\) to be numerical.

## 2. Exact spectral and Bochner ledger

Isotropy and the weak eigenvalue equation give

\[
\begin{aligned}
 &\mathbb Er=0,\qquad \mathbb E[Xr]=0,\qquad
 \|r\|_2^2=\delta^2,\\
 &\mathbb E\nabla f=\lambda a,\qquad
 \mathbb E\nabla r=-(1-\lambda)a,\\
 &\mathbb E|\nabla r|^2
   =1-\lambda+(2\lambda-1)\delta^2 .
\end{aligned}                                                \tag{2.1}
\]

For smooth full-support measures, with the standard Reilly boundary term
for convex supports, put

\[
 B_f=\mathbb E\|D^2f\|_{\rm HS}^2,\qquad
 \mathcal R_V(f)\ge0
\]

for the sum of the potential-curvature and boundary-curvature terms.
Bochner--Reilly and componentwise Poincaré yield the exact identities

\[
 B_f+\mathcal R_V(f)=\lambda^2,                       \tag{2.2}
\]

\[
 \operatorname {Var}(\nabla f)
 =\lambda\big[1-\lambda(1-\delta^2)\big],             \tag{2.3}
\]

and the nonnegative defect decomposition

\[
 \boxed{\;
 [B_f-\lambda\operatorname {Var}(\nabla f)]
 +\mathcal R_V(f)
 =\lambda^3(1-\delta^2).\;}                           \tag{2.4}
\]

In particular \(\|D^2r\|_2=\|D^2f\|_2\le\lambda\).  These formulas contain
no unproved KLS input.

The tempting interpolation

\[
 \operatorname {Var}(\nabla f)
 \le C\,\delta\,\|D^2f\|_2
\]

is false.  The first degree-one Neumann mode on the isotropic Euclidean
ball satisfies

\[
 \lambda=1-\frac1n+O(n^{-2}),\quad
 \delta=\frac1n+O(n^{-2}),\quad
 B_f=\frac3n+O(n^{-2}),\quad
 \operatorname {Var}(\nabla f)=\frac1n+O(n^{-2}),
\]

so the ratio grows like \(\sqrt{n/3}\).

## 3. The surviving mean-gradient gate

The narrower statement is

\[
 \boxed{\quad
 |\mathbb E\nabla g|
 \le C\big(\|g\|_2+\|D^2g\|_2\big),
 \qquad
 \mathbb Eg=0,\quad \mathbb E[Xg]=0 .
 \quad}                                                \tag{MG}
\]

Applied to \(r\), it gives

\[
 (1-\lambda)\sqrt{1-\delta^2}
 \le C(\delta+\lambda),                               \tag{3.1}
\]

which forces \(\lambda\ge c(C)>0\) once \(\delta\) is a sufficiently small
universal number.

There is a complete one-dimensional proof of the stronger Hessian-only
form.  If \(\nu\) is centered, variance one, and log-concave, its
one-dimensional Poincaré constant is at most \(12\).  For
\(h\perp\{1,x\}\), setting \(m=\int h'\,d\nu\) and applying Poincaré to
\(h'\) and then to \(h-mx\) gives

\[
 \left|\int h'\,d\nu\right|
 \le12\left(\int(h'')^2d\nu\right)^{1/2}.             \tag{3.2}
\]

For a product of such laws, the conditional projections
\(g_i(x_i)=\mathbb E[g\mid X_i=x_i]\) are orthogonal and

\[
\begin{aligned}
 |\mathbb E\nabla g|^2
 &\le144\sum_i\|g_i''\|_2^2\\
 &\le144\,\mathbb E\|D^2g\|_{\rm HS}^2 .
\end{aligned}                                        \tag{3.3}
\]

Thus there is no hidden dimension loss for products.

For a symmetric matrix \(B\), let

\[
 q_B=X^TBX-\operatorname {tr}B,\qquad
 T B=\mathbb E[Xq_B],\qquad
 g_B=q_B-(TB)\cdot X .
\]

If \(\mathsf C\) is the covariance operator of centered quadratic forms and
\(\mathsf S=\mathsf C-T^*T\), then

\[
 \|g_B\|_2^2=\langle B,\mathsf S B\rangle,\qquad
 \mathbb E\nabla g_B=-TB,\qquad D^2g_B=2B.            \tag{3.4}
\]

The Hessian-only quadratic gate is exactly

\[
 \sup_{|u|=1}
 \left\|\mathbb E[(u\cdot X)(XX^T-I)]\right\|_{\rm HS}
 \le C.                                               \tag{3.5}
\]

Projection thin shell gives, for every rank-\(k\) projection \(P\),

\[
 |\operatorname {tr}(PM_u)|\le C\sqrt{k}.
\]

Hence \(|\lambda_j(M_u)|\le Cj^{-1/2}\) and only
\(\|M_u\|_{\rm HS}\le C\sqrt{\log(en)}\).  The harmonic diagonal matrix
shows that these rank-by-rank data alone cannot remove the logarithm.

Disintegrating along \(u\) also does not prove (MG).  For
\(h(t)=\mathbb E[g\mid u\cdot X=t]\) and conditional score
\(s_t=V_t-\mathbb E_tV_t\),

\[
\begin{aligned}
 h'&=\mathbb E_tg_t-\mathbb E_t(gs_t),\\
 h''&=\mathbb E_tg_{tt}-2\mathbb E_t(g_ts_t)
 +\mathbb E_t[g(s_t^2-\partial_ts_t)] .
\end{aligned}                                        \tag{3.6}
\]

Prékopa controls only the scalar Fisher difference in (3.6); it does not
control the signed score covariances.  For hard supports these terms are
moving-boundary fluxes.

## 4. Dependent classes where the gate closes

For the paired box wedge

\[
 K_a=\{(t,x,y):|t|\le1,\ |x_i|\le1+a_it,\
                         |y_i|\le1-a_it\},
\]

put \(S=\sum_i a_i^2\) and \(v=\operatorname {Var}T\).  The axial
slice potential satisfies

\[
 W''(t)\ge2S,
\qquad\text{hence}\qquad vS\le\frac12 .               \tag{4.1}
\]

After exact covariance whitening, the only third coefficients are

\[
 \lambda_i=\frac{2a_i\sqrt v}{1+a_i^2v},
\qquad \sum_i\lambda_i^2\le2,
\]

and therefore

\[
 \sup_{|u|=1}\|M_u\|_{\rm HS}\le2.                    \tag{4.2}
\]

The whitened law is the image of a product of isotropic one-dimensional
log-concave laws under an explicit map \(F\) with
\(\|DF\|_{\rm op}\le4\).  Consequently

\[
 C_P(K_a^{\rm iso})\le192,\qquad
 |\mathbb E\nabla g|\le192\|D^2g\|_2
\quad(g\perp\mathrm {Aff}).                           \tag{4.3}
\]

A genuinely one-sided wedge has

\[
 \|DF\|_{\rm op}\le\sqrt{28}+3,\qquad
 C_P\le444+72\sqrt{28}<972,                           \tag{4.4}
\]

and obeys the same type of Hessian-only estimate.  Gaussian convolution
preserves these dimension-free conclusions and produces smooth,
full-support, genuinely dependent examples.

Several further exact classes close the quadratic gate:

* For every log-concave Dirichlet law with parameters
  \(\alpha_i\ge1\), \(A=\sum_i\alpha_i\), put
  \(q_i=\alpha_i/A\).  A unit tangent direction is represented by
  coefficients \(a_i\) satisfying
  \(\sum_iq_ia_i=0\) and \(\sum_iq_ia_i^2=1\), and
  \[
  \|M_{u(a)}\|_{\rm HS}^2
  =\frac{4(A+1)}{(A+2)^2}
    \left(\sum_i a_i^2-2\right)<4.                    \tag{4.5}
  \]
  The constant \(2\) is asymptotically sharp.
* For cones over centrally symmetric bases,
  \[
  \sup_{|u|=1}\|M_u\|_{\rm HS}^2
  =\frac{4(n-1)(n+2)}{(n+3)^2}\longrightarrow4.       \tag{4.6}
  \]
* Affine box slices and affine positive-definite Gaussian covariance
  pencils have universal bounds.  In both cases the second derivative of
  the axial log-volume charges precisely the squared Frobenius norm of the
  normalized covariance velocity.

These results rule out the most direct realization of the harmonic
\(j^{-1/2}\) spectrum.  What remains uncontrolled is simultaneous
non-affine shape change and noncommuting rotation of conditional
covariances.

## 5. Function-specific localization and convex-set diagnostics

For a fixed normalized first eigenfunction, a hard localization driver
with \(C_tb_t=0\), \(b_t=\operatorname {Cov}_t(f,X)\), preserves its
posterior mean pathwise.  If

\[
 \mathcal E_t=\mathbb E_t|\nabla f|^2,\qquad
 R_T=\text{terminal survivor variance},
\]

the rank-one defect theorem gives the exact ledger

\[
 1\le\frac{192}{T}\lambda
      +96\,\mathbb E[R_T\mathcal E_T].                \tag{5.1}
\]

Under the energy-biased law
\(dQ=(\mathcal E_T/\lambda)dP\), a small gap therefore forces

\[
 \mathbb E_QR_T\ge\frac1{96\lambda}-\frac2T .         \tag{5.2}
\]

The soft driver creates full curvature \(\delta TI\) only at the
complementary cost

\[
 1-\mathbb Em_T^2\le\frac{\lambda}{\delta T}.         \tag{5.3}
\]

Thus this scheme cannot preserve the signal and regularize its selected
direction at a universal cost without a new adaptive-survivor theorem.

For a closed convex half-set \(A\), \(\mu(A)=1/2\), put
\(p=\mu^+(A)\).  Log-concavity of
\(F(t)=\mu(A+tB_2^n)\) gives the sharp scalar implication

\[
 \boxed{\quad
 \mathbb E\,d(X,A)\ge
 \frac{2\log2-1}{4p}.
 \quad}                                                \tag{5.4}
\]

Conversely, the distance function shows
\(C_P(\mu)\ge(\mathbb E d(X,A))^2\).  A universal upper bound in (5.4)
is therefore itself KLS-strength.  Thin shell handles centered balls and
isotropy handles halfspaces, but neither controls arbitrary convex
half-sets.  Naive convexification of a general Cheeger witness fails even
for the isotropic Laplace law.

## 6. Displacement audit and next round

For every balanced Brenier map, the midpoint Cayley potential \(q\) is
globally \(1\)-smooth and satisfies

\[
 \frac18W_2^2\le\mathbb E_\mu[\sigma q]\le\frac38W_2^2,
\qquad
 \int|\nabla q|^2d\mu\le W_2^2.                       \tag{6.1}
\]

A Poincaré bound restricted to all such generally nonconvex Cayley
potentials is quantitatively equivalent to KLS.  In the firm zero-strain
subbranch \(q\) is convex.  However, the outer-halves-to-middle-halves
interval transport has zero entropy deficit and

\[
 \mathbb E_{\rm source}r\ge\mathbb E_{\rm target}r
\]

for every convex \(r\).  Hence no prescribed-orientation signed convex
certificate, even one chosen from the transport, can recover the positive
quadratic displacement with an entropy-deficit error.  The nonconvex
Cayley potential is essential.

The next independent routes are:

1. exploit the special eigenfunction equation to prove (MG) without
   proving the full directional third-moment tensor bound;
2. extend slice-curvature compensation from commuting affine pencils to
   noncommuting conditional shape flows;
3. test whether the parallel coupling of log-affine tilts controls a low
   spectral mode rather than only the trace
   \(\sum_i\|x_i\|_{H^{-1}}^2\);
4. seek a global laminar or incidence theorem for the cyclic
   zero-strain displacement branch.

Every item remains independent of the ordinary
\(\sqrt{\log n}\) covariance bootstrap.  None is presently a completed
KLS proof.
