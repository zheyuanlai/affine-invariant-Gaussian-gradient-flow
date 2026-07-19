# Clean-room audit of the smooth Gaussian-fiber centroid and WFI bounds

## 0. Verdict

The dimension-free centroid theorem in Section B of
shape_centroid_stability_audited.md is correct for the smooth Gaussian-fiber
class, after one repair to its one-dimensional input.  The displayed
pointwise assertion

\[
        \tau(s)|\Phi'(s)|\leq 1+|s|
\]

does not follow from the tail calculation given there: if
\(a=|\Phi'(s)|\), that calculation gives
\(|s|+1/a\), not \(|s|+1\).  No counterexample to the assertion itself is
needed.  The theorem only needs a universal constant in front, and the
following clean-room proof establishes

\[
        \tau(s)|\Phi'(s)|\leq 400(1+|s|).
        \tag{0.1}
\]

With this replacement, every subsequent constant remains universal and
independent of the fiber dimension.

The full Gaussian-fiber weighted-Fisher calculation in
gaussian_full_wfi.md is also correct in the smooth setting.  With the
generator convention stated below, the exact Sylvester equation is

\[
        QA+AQ=Q'=-QR'Q,
        \qquad\text{equivalently}\qquad
        AR+RA=-R'.
\]

The sharp comparisons used there are

\[
 C_s\leq \frac12J,\qquad B_s^0\leq J,\qquad J\leq\Phi''.
\]

Consequently the shape part costs at most \(3/2\), and the assumption
\(R(s)\succeq I\) converts the centroid curvature term into the corrected
centroid estimate.  The claimed nonsmooth/distributional extension is not
proved by either note and must remain a separate approximation problem.

## 1. Exact hypotheses

Let \(J\subseteq\mathbb R\) be an interval and let

\[
 U(s,z)=W(s)+\frac12(z-m(s))^TQ(s)(z-m(s)),
 \qquad R(s)=Q(s)^{-1}\succ0,
 \tag{1.1}
\]

where \(W:J\to\mathbb R\), \(m:J\to\mathbb R^d\), and
\(R:J\to\mathbb S_{++}^d\) are \(C^2\) in the interior of \(J\).
Assume the extended-valued potential, equal to \(+\infty\) off
\(J\times\mathbb R^d\), is convex and \(e^{-U}\) is normalizable.  Its
\(S\)-marginal is

\[
 \rho(s)=e^{-\Phi(s)},\qquad
 \Phi(s)=W(s)-\frac12\log\det R(s)+\text{constant}.
 \tag{1.2}
\]

Assume

\[
 \mathbb ES=0,\qquad \mathbb ES^2=1,\qquad
 \mathbb E[S\,m(S)]=0,\qquad
 \mathbb ER(S)\preceq\Lambda I_d.
 \tag{1.3}
\]

The vector expectation in (1.3) is assumed absolutely integrable.
For the full WFI conclusion, assume in addition

\[
                         R(s)\succeq I_d
                         \quad\text{for all }s\in J.
 \tag{1.4}
\]

The proof below is for these smooth hypotheses.  Hard endpoints are allowed
provided derivatives are interpreted one-sidedly.  A limit of smooth
families is covered only if convergence of every displayed nonnegative
quadratic form and preservation of (1.3)-(1.4) are separately verified.

## 2. One-dimensional lemmas

Let \(\rho=e^{-\Phi}\) be any centered variance-one log-concave density on
an interval, and let

\[
 \tau(s)\rho(s)=\int_s^{\sup J}t\rho(t)\,dt
              =-\int_{\inf J}^s t\rho(t)\,dt
 \tag{2.1}
\]

be its canonical Stein kernel.

### 2.1 Density, mode, and moment bounds

The elementary one-dimensional isotropic density bounds give

\[
 \rho(0)\geq\frac1{20},\qquad \|\rho\|_\infty\leq1.
 \tag{2.2}
\]

For completeness, the lower bound follows from Grünbaum's
\(\mathbb P(S\geq0),\mathbb P(S\leq0)\geq e^{-1}\), Chebyshev's
\(\mathbb P(|S|\geq2)\leq1/4\), and unimodality.  On the side of zero
opposite a mode, \(\rho\leq\rho(0)\), while that side inside \([-2,2]\)
has mass at least \(e^{-1}-1/4\).  Hence
\(2\rho(0)\geq e^{-1}-1/4>1/10\).  The upper bound in (2.2) is the standard
one-dimensional log-concave extremal inequality
\(\|\rho\|_\infty^2\operatorname {Var}(S)\leq1\), whose extremal
one-sided law is exponential.

If \(s_0\) is a mode, log-concavity gives
\(\rho(t)\geq\rho(0)\) on the segment joining \(0\) to \(s_0\).
Thus

\[
                             |s_0|\leq20.
 \tag{2.3}
\]

Borell's one-dimensional moment lemma gives a numerical \(M_3\) such that

\[
                             \mathbb E|S|^3\leq M_3.
 \tag{2.4}
\]

Only finiteness of this universal number is used below.

### 2.2 Concave core/growth lemma

If \(h:J\to[0,\infty)\) is concave and
\(\mathbb Eh(S)\leq\Lambda\), then

\[
                         h(s)\leq1000\Lambda(1+|s|).
 \tag{2.5}
\]

Here is a proof not using a pointwise core assertion.  Let \(q_\alpha\)
denote the \(\alpha\)-quantile.  Markov's inequality shows that
\(\{h\leq20\Lambda\}\) has probability at least \(19/20\).  Hence each of

\[
 [q_{.1},q_{.2}],\quad[q_{.3},q_{.4}],\quad
 [q_{.6},q_{.7}],\quad[q_{.8},q_{.9}]
\]

contains a point at which \(h\leq20\Lambda\).  By
\(\|\rho\|_\infty\leq1\), points chosen in the first two intervals are
separated by at least \(0.1\); the same holds for points in the last two
intervals.  Chebyshev places all four points in
\([-\sqrt {10},\sqrt {10}]\).  Monotonicity of secant slopes of a concave
function and nonnegativity of \(h\) then give, to the right of the first
pair and to the left of the second pair,

\[
 h(s)\leq20\Lambda+200\Lambda(|s|+\sqrt {10}).
\]

These two regions cover \(J\), and (2.5) follows.

### 2.3 Correct Stein-kernel and score bounds

Let \(T_+(s)=\mathbb P(S\geq s)\).  Log-concavity of \(\rho\) implies
log-concavity of \(T_+\), so its hazard
\(\lambda_+(s)=\rho(s)/T_+(s)\) is nondecreasing.  By (2.2), for \(s\geq0\),

\[
 {T_+(s)\over\rho(s)}
 \leq {T_+(0)\over\rho(0)}\leq20.
 \tag{2.6}
\]

The tangent bound for the log-concave survival function gives

\[
 \int_s^{\sup J}T_+(t)\,dt
 \leq {T_+(s)\over\lambda_+(s)}
 ={T_+(s)^2\over\rho(s)}.
 \tag{2.7}
\]

Integration by parts in (2.1), followed by (2.6)-(2.7), yields

\[
 \tau(s)
 ={sT_+(s)+\int_s^{\sup J}T_+(t)\,dt\over\rho(s)}
 \leq20s+400\qquad(s\geq0).
\]

The lower-tail argument is identical.  Therefore

\[
                   0\leq\tau(s)\leq400(1+|s|),
 \qquad
                   \mathbb E\tau(S)^2\leq640000.
 \tag{2.8}
\]

Now put \(a=|\Phi'(s)|\) at a differentiability point.  If \(a\geq1\)
and \(\Phi'(s)>0\), convexity gives
\(\rho(s+u)\leq\rho(s)e^{-au}\), and the right-tail formula gives

\[
 \begin{aligned}
 a\tau(s)
 &\leq {a\over\rho(s)}
       \int_0^\infty (|s|+u)\rho(s+u)\,du\\
 &\leq |s|+{1\over a}\leq1+|s|.
 \end{aligned}
\]

The left-tail case \(\Phi'(s)<0\) is the same.  If \(a<1\), use (2.8).
Thus the corrected universal bound (0.1) holds.  In particular,

\[
 \mathbb E\!\left[
   \tau(S)^2(1+|S|)\Phi'(S)^2\right]
 \leq 400^2\,\mathbb E(1+|S|)^3
 \leq 4\cdot400^2(1+M_3).
 \tag{2.9}
\]

The constant \(1\) in the original pointwise assertion is not obtained by
its proof; (0.1) is the rigorously established replacement needed below.

### 2.4 Stein-curvature identity

At interior differentiability points, (2.1) gives

\[
                         \tau'-\tau\Phi'=-s.
 \tag{2.10}
\]

The bounds above justify integration by parts with
\(g=\tau^2\Phi'\).  Indeed, \(\tau\Phi'\) and \(\tau'\) have at most
linear growth, while one-dimensional isotropic log-concave laws have
uniform exponential tails; the boundary term
\(\rho\tau^2\Phi'\) vanishes at both endpoints.  Hence

\[
 \begin{aligned}
 \mathbb E[\tau^2\Phi'']
 &=\mathbb E[\tau^2(\Phi')^2-2\tau\tau'\Phi']\\
 &=\mathbb E[S^2-(\tau')^2]
 =1-\mathbb E(\tau')^2\leq1.
 \end{aligned}
 \tag{2.11}
\]

For a nonsmooth convex \(\Phi\), the same statement uses its curvature
measure and requires a separate monotone-approximation argument.  Only the
smooth identity (2.11) is used here.

## 3. Clean-room proof of the Gaussian centroid theorem

We first compute the exact joint-convexity constraint.  Put
\(y=z-m(s)\).  At fixed \(z\),

\[
 \begin{aligned}
 U_{zz}&=Q,\\
 U_{sz}&=Q'y-Qm',\\
 U_{ss}&=W''+\frac12y^TQ''y-2m'^TQ'y
          -m''{}^TQy+m'^TQm'.
 \end{aligned}
 \tag{3.1}
\]

The Schur complement of \(Q\) is

\[
 W''-m''{}^TQy+
 \frac12y^T\!\left(Q''-2Q'R Q'\right)y\geq0.
\]

Twice differentiating \(R=Q^{-1}\) gives

\[
                   Q''-2Q'RQ'=Q(-R'')Q.
\]

Thus, with \(K=-R''\) and \(x=Qy\),

\[
                   W''-m''\!\cdot x+\frac12x^TKx\geq0
                   \quad\text{for every }x\in\mathbb R^d.
 \tag{3.2}
\]

The quadratic polynomial in (3.2) is nonnegative for every \(x\).
Consequently

\[
 K\succeq0,\qquad m''\in\operatorname {Ran}K,\qquad
 m''{}^TK^\dagger m''\leq2W''.
 \tag{3.3}
\]

For every unit vector \(u\), Cauchy-Schwarz in the \(K\)-seminorm gives

\[
 |u\cdot m''|^2
 \leq2W''\,u^TKu
 \leq2W''\,\lambda_{\max}(R)\operatorname {tr}(R^{-1}K).
 \tag{3.4}
\]

The last inequality follows from
\(R^{-1}\succeq\lambda_{\max}(R)^{-1}I\) and
\(u^TKu\leq\operatorname {tr}K\).

Differentiating (1.2) gives the exact marginal-curvature formula

\[
 \Phi''
 =W''+\frac12\operatorname {tr}(R^{-1}K)
 +\frac12
  \|R^{-1/2}R'R^{-1/2}\|_{\mathrm {HS}}^2.
 \tag{3.5}
\]

All three terms on the right are nonnegative.  Since
\(2\sqrt{xy}\leq x+y\), (3.4)-(3.5) imply the slightly sharper version

\[
                         |m''(s)|
 \leq\sqrt{\lambda_{\max}R(s)}\,\Phi''(s).
 \tag{3.6}
\]

For each unit \(u\), the scalar function \(u^TR(s)u\) is nonnegative and
concave by \(K\succeq0\), and its expectation is at most \(\Lambda\).
Applying (2.5) and then taking the supremum over \(u\) gives

\[
                 \lambda_{\max}R(s)
 \leq1000\Lambda(1+|s|).
 \tag{3.7}
\]

Choose a mode \(s_0\) of \(\rho\), and let \(b\) denote the appropriate
one-sided value of \(m'(s_0)\).  If the mode is interior,
\(\Phi'(s_0)=0\); at a hard endpoint, convexity gives the same inequality
below because the interior curvature accumulated from the endpoint is at
most \(|\Phi'(s)|\).  By (2.3), \(|s_0|\leq20\).  Integrating (3.6) on
the segment from \(s_0\) to \(s\), using monotonicity of \(\Phi'\), gives

\[
 |m'(s)-b|
 \leq\sqrt{21000\,\Lambda}\,
      \sqrt{1+|s|}\,|\Phi'(s)|.
 \tag{3.8}
\]

Combining (3.8) with (2.9),

\[
 \begin{aligned}
 D&:=\mathbb E[\tau(S)^2|m'(S)-b|^2]\\
 &\leq
 21000\cdot4\cdot400^2(1+M_3)\,\Lambda
 =:C_2\Lambda.
 \end{aligned}
 \tag{3.9}
\]

The canonical Stein identity is applicable componentwise because (3.8)
and (2.8) give the needed integrability.  From (1.3),

\[
              0=\mathbb E[S\,m(S)]
               =\mathbb E[\tau(S)m'(S)].
 \tag{3.10}
\]

Also \(\mathbb E\tau=\mathbb ES^2=1\).  Therefore

\[
 b=-\mathbb E[\tau(m'-b)],\qquad |b|^2\leq D.
 \tag{3.11}
\]

Finally, (2.8), (3.9), and (3.11) give

\[
 \boxed{\;
 \mathbb E[\tau(S)^2|m'(S)|^2]
 \leq2C_2(1+640000)\Lambda
 =:C_{\rm cent}\Lambda .
 \;}
 \tag{3.12}
\]

This proves the claimed dimension-free Gaussian centroid theorem, with a
deliberately loose but explicit universal coefficient.  No conditional
Poincare inequality and no lower bound on \(R\) entered this part.

## 4. Clean-room proof of the full Gaussian-fiber WFI theorem

Assume now (1.4).  For fixed \(s\), let

\[
 L_s=\Delta_z-\langle Q(s)(z-m(s)),\nabla_z\rangle
\]

be the nonpositive Ornstein-Uhlenbeck generator on
\(N(m(s),R(s))\), and put \(\ell_s=\partial_s\log q_s\).  Direct
differentiation gives

\[
 \ell_s
 =m'^TQy-\frac12
   \left(y^TQ'y-\operatorname {tr}(Q'R)\right).
 \tag{4.1}
\]

There is a unique symmetric solution \(A=A(s)\) of

\[
 QA+AQ=Q'=-QR'Q,
 \qquad\text{equivalently}\qquad
 AR+RA=-R'.
 \tag{4.2}
\]

The function

\[
                 g_s(z)=-m'(s)^Tz+\frac12y^TAy
\]

solves \(L_sg_s=\ell_s\), up to an irrelevant additive constant.  Indeed,
\(\operatorname {tr}A=\frac12\operatorname {tr}(Q'R)\), which follows by
multiplying (4.2) by \(R\) and taking the trace.  Hence

\[
 F_s:=\nabla_zg_s=-m'+Ay,\qquad
 C_s:=\mathbb E_s\|D_zF_s\|_{\rm HS}^2=\|A\|_{\rm HS}^2,
 \tag{4.3}
\]

and

\[
 \begin{aligned}
 B_s&:=\mathbb E_s\langle QF_s,F_s\rangle\\
 &=m'^TQm'+B_s^0,\qquad
 B_s^0=\operatorname {tr}(QARA).
 \end{aligned}
 \tag{4.4}
\]

Diagonalize \(R\) at this fixed value of \(s\); no derivative of the
eigenbasis is taken.  If its eigenvalues are \(r_i>0\), then

\[
 A_{ij}=-{R'_{ij}\over r_i+r_j}.
 \tag{4.5}
\]

Define

\[
 J_s:=\frac12\operatorname {tr}
       (R^{-1}R'R^{-1}R')
     =\frac12\sum_{i,j}{(R'_{ij})^2\over r_ir_j}.
 \tag{4.6}
\]

For diagonal terms, the contribution to \(B_s^0\) is one half of the
corresponding contribution to \(J_s\).  For \(i<j\), the two ordered
terms in \(B_s^0\) combine to

\[
 (R'_{ij})^2{r_i^2+r_j^2\over
                  r_ir_j(r_i+r_j)^2}
 \leq{(R'_{ij})^2\over r_ir_j},
\]

which is the pair contribution to \(J_s\).  Thus

\[
                              0\leq B_s^0\leq J_s.
 \tag{4.7}
\]

Likewise, the diagonal ratio for \(C_s\) is \(1/2\), while for an
off-diagonal pair it is

\[
 {2r_ir_j\over(r_i+r_j)^2}\leq\frac12.
\]

Therefore

\[
                              0\leq C_s\leq\frac12J_s.
 \tag{4.8}
\]

Formula (3.5) says exactly

\[
 \Phi''=W''+\frac12\operatorname {tr}(R^{-1}(-R''))+J_s.
\]

The Schur-complement conclusion (3.3) makes the first two terms
nonnegative, so

\[
                             J_s\leq\Phi''(s).
 \tag{4.9}
\]

Equations (4.7)-(4.9) and the Stein-curvature identity (2.11) give

\[
 \int_J\rho\tau^2(C_s+B_s^0)\,ds
 \leq\frac32\int_J\rho\tau^2\Phi''\,ds
 \leq\frac32.
 \tag{4.10}
\]

Finally, \(R\succeq I\) implies \(Q\preceq I\).  The centroid theorem
(3.12) therefore gives

\[
 \int_J\rho\tau^2m'^TQm'\,ds
 \leq\mathbb E[\tau^2|m'|^2]
 \leq C_{\rm cent}\Lambda.
 \tag{4.11}
\]

Combining (4.4), (4.10), and (4.11) proves

\[
 \boxed{\;
 \int_J\rho(s)\tau(s)^2(C_s+B_s)\,ds
 \leq\frac32+C_{\rm cent}\Lambda .
 \;}
 \tag{4.12}
\]

Every matrix comparison is pointwise and independent of \(d\); no
commutativity of \(R\) and \(R'\) is used.

## 5. Scope and gaps that must remain explicit

1. The exact constant-one pointwise score inequality in the audited note is
   not proved by its tail argument.  The corrected constant-400 estimate is
   sufficient for the theorem.

2. The centroid proof requires the mixed moment in (1.3) to exist
   absolutely so that the Stein cancellation is meaningful.

3. Equation (2.11) is proved here for a smooth marginal with the stated
   boundary decay, which follows from the one-dimensional tail estimates.
   Curvature atoms and hard nonsmooth limits require a monotone
   approximation argument.

4. The lower covariance bound \(R\succeq I\) is essential for the last
   comparison used in this proof.  More generally \(R\succeq\kappa I\)
   gives the centroid cost \(C_{\rm cent}\Lambda/\kappa\).  An averaged
   lower bound does not suffice for this step.

5. The assertion that a distributional \(R\) is covered whenever it is
   "obtained by smooth approximation" is not itself an approximation
   theorem.  One must prove preservation of joint convexity, the pointwise
   lower covariance bound, the expected upper covariance bound, and
   lower-semicontinuity or convergence of \(C_s\) and \(B_s\).

6. Nothing here extends the theorem to non-Gaussian conditional fibers or
   removes unit transverse Gaussian noise.  Those are the remaining
   conjecture-strength steps.
