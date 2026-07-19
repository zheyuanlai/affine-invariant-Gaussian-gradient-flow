# The heat-generated physical selector: analytic replacement, information bounds, and a phase no-go

## 0. Outcome

Let the fixed-scale construction of
fixed_scale_physical_splicing.md be in force.  Thus

\[
 \alpha=10^{-10},\qquad s=\alpha K,\qquad
 \Delta:=\mathbb E(R-|m|)<\epsilon_*p,\qquad
 \epsilon_*=6.02\cdot10^{-5},                         \tag{0.1}
\]

and the good posterior flux satisfies

\[
 b:=\mathbb E r_G\ge b_0p,\qquad
 b_0={.13(1-10^{-5})\over8\pi}>.00517.               \tag{0.2}
\]

Here

\[
 R(x)=P_s|W|(x),\qquad r_G(x)=P_s(1_G|W|)(x),\qquad
 m(x)=P_sW(x)=\nabla F_0(x).                         \tag{0.3}
\]

The first result is a strict improvement of the hard selector used in the
fixed-scale report.  Define

\[
 q(x)={|m(x)|\over R(x)},\qquad
 \omega_{\rm an}(x)={r_G(x)\over R(x)},\qquad
 \eta(x)=\omega_{\rm an}(x)|m(x)|=r_G(x)q(x).        \tag{0.4}
\]

The scalar selector \(\omega_{\rm an}\) is analytic, belongs to \([0,1]\),
and gives a genuine physical coarea submeasure.  Its normal matrix satisfies

\[
 \boxed{\operatorname{tr}M_{\rm an}>.005109p,\qquad
        {\operatorname{tr}M_{\rm an}\over
         \|M_{\rm an}\|_{op}}>19.54.}                \tag{0.5}
\]

Moreover, restriction to

\[
 \mathcal C=\{\omega_{\rm an}\ge1/2000,\ q\ge1/2\}   \tag{0.6}
\]

may be done in two stages.  Restricting only to
\(\{\omega_{\rm an}\ge1/2000\}\) leaves trace \(>.004609p\) and effective
rank \(>17.62\).  The smaller core \(\mathcal C\) still leaves

\[
 \boxed{\operatorname{tr}M_{\mathcal C}>.004489p,\qquad
        {\operatorname{tr}M_{\mathcal C}\over
         \|M_{\mathcal C}\|_{op}}>17.16.}            \tag{0.7}
\]

Thus the hard set \(H=\{d\le\lambda r_G\}\) is not intrinsic.  On the
retained physical flux, a fixed fraction of the \(|W|\)-tilted Gaussian
channel is good, and the total directional mean has magnitude at least
one half.

The second result is a collection of dimension-free Gaussian-scale
variation estimates.  If

\[
 B_0=\mathbb E R=\mathbb E|W|\le p,\qquad
 M_s={I_0\over\sqrt s},\qquad I_0={1\over\sqrt{2\pi}},             \tag{0.8}
\]

then

\[
 \int R\left\|\nabla{m\over R}\right\|_{HS}^2d\mu
 \le {8B_0\over s}\left(1+\log{M_s\over B_0}\right),              \tag{0.9}
\]

\[
 \int r_G\|\nabla v_G\|_{HS}^2d\mu
 \le {8b\over s}\left(1+\log{M_s\over b}\right),\qquad
 v_G={P_s(1_GW)\over r_G},                                      \tag{0.10}
\]

and the analytic selector has the binary Fisher bound

\[
 \boxed{
 \int R\,{|\nabla\omega_{\rm an}|^2
                 \over\omega_{\rm an}(1-\omega_{\rm an})}\,d\mu
 \le {2B_0\over s}
       \left(\log{M_s\over B_0}+\log2\right).}        \tag{0.11}
\]

The usual zero conventions apply at \(\omega=0,1\).  By Buser--Ledoux and
the half-mass characterization of the Cheeger constant,
\(p\sqrt K\ge c_{\rm BL}/2\).  Hence every logarithm in (0.9)--(0.11) is at
most

\[
 L_\alpha=\log {2I_0\over b_0c_{\rm BL}\sqrt\alpha},              \tag{0.12}
\]

a fixed universal number once \(\alpha=10^{-10}\) is frozen.

The third result transfers a between-level normal-projector variance to the
smooth good-direction field.  If the between-level term of the analytic
physical measure is at least \(8/17\), then

\[
 \boxed{
 \operatorname{Var}\!\left(
   \mathbb E_\nu[v_Gv_G^T\mid F_0]\right)>0.143,}     \tag{0.13}
\]

where \(d\nu=\eta\,d\mu/\int\eta d\mu\).  Thus the heat origin does prevent
the completely arbitrary rotating-cap selector: a fixed between-level
charge must also occur in an explicitly smoothed conditional direction.

It still does not yield a \(10^{-4}p\) Bernstein or coarea charge.  The
available estimates are upper energy bounds.  Turning (0.13) into a lower
Dirichlet bound requires a capacity inequality for the selected physical
measure.  Even optimistically using the ambient Poincare scale \(K=s/\alpha\)
only gives energy of order

\[
                         {p\over K}={\alpha p\over s}.             \tag{0.14}
\]

After multiplication by the natural Gaussian factor \(s\), this is
\(\alpha p=10^{-10}p\), a factor \(10^6\) below the required
\(10^{-4}p\).  This is exactly the old scalar power obstruction, now exposed
after the selector has been fully regularized.

Finally, Section 7 gives an exact local heat-generated phase profile in a
log-affine chart.  Several flat facets with different normal slopes produce
a softmax of normal projectors across physical levels.  Each facet has
\(e=1\), affine profile coordinate, zero Hessian charge, zero eikonal charge,
and the analytic selector (0.4).  The only possible charge is at the ridges,
ends, or global transition region.  Therefore Gaussian derivative and
information bounds alone cannot prove the desired splice lemma; a physical
incidence or an extremality-dependent capacity theorem is still necessary.

## 1. Setup and the analytic physical selector

Let \(X\sim\mu\), let \(Z\sim\gamma_n\) be independent, and write

\[
                         P_sh(x)=\mathbb E h(x+\sqrt sZ).           \tag{1.1}
\]

At the fixed heat time, let

\[
 W(y)=\nabla g(y),\qquad
 u(y)={W(y)\over|W(y)|},\qquad
 h_G(y)=1_G(y)|W(y)|.                              \tag{1.2}
\]

The Gaussian centroid inequality gives

\[
                         |W|\le M_s={I_0\over\sqrt s}.             \tag{1.3}
\]

Define \(R,r_G,m,q,\omega_{\rm an},\eta\) by (0.3)--(0.4).  Since
\(0\le h_G\le |W|\),

\[
 0<r_G\le R,\qquad 0\le q\le1,\qquad
 0<\omega_{\rm an}\le1.                             \tag{1.4}
\]

Strict positivity follows from \(b>0\) and positivity of the Gaussian
kernel.  Both \(R\) and \(r_G\) are real analytic even though \(G\) is
defined by hard inequalities.  Thus \(\omega_{\rm an}=r_G/R\) is real
analytic.

Let

\[
 \theta={m\over|m|}
\]

where \(m\ne0\), with an arbitrary value on the zero set.  The analytic
physical matrix is

\[
 M_{\rm an}=\int\eta(x)\theta(x)\theta(x)^T\,d\mu(x).              \tag{1.5}
\]

Because \(\eta=\omega_{\rm an}|m|\) and
\(0\le\omega_{\rm an}\le1\), weighted coarea gives

\[
 M_{\rm an}=\int_0^1\int_{\partial^*A_r}
 \omega_{\rm an}\,n_rn_r^Te^{-V}\,
 d\mathcal H^{n-1}\,dr.                             \tag{1.6}
\]

Thus it is a genuine submeasure of the physical level boundaries.

### 1.1 Trace and rank

Recall the matrix

\[
 D=\int r_G\theta\theta^T\,d\mu,\qquad \operatorname{tr}D=b.       \tag{1.7}
\]

The exact alignment calculation in the fixed-scale report gives

\[
                         \|D\|_{op}\le2\|B_G\|_{op}+4\Delta,       \tag{1.8}
\]

where \(B_G\) is the good posterior flux matrix and

\[
                         {b\over\|B_G\|_{op}}\ge500.               \tag{1.9}
\]

Since \(R-|m|=R(1-q)=d\),

\[
 r_G-\eta=r_G(1-q)={r_G\over R}d\le d.              \tag{1.10}
\]

Consequently

\[
 0\preceq D-M_{\rm an},\qquad
 \operatorname{tr}(D-M_{\rm an})\le\Delta.          \tag{1.11}
\]

It follows that

\[
 \operatorname{tr}M_{\rm an}\ge b-\Delta,\qquad
 \|M_{\rm an}\|_{op}\le {2b\over500}+4\Delta.        \tag{1.12}
\]

Both ratios in (1.12) are monotone in \(b\) in the direction needed below.
Using \(b\ge b_0p\) and \(\Delta<\epsilon_*p\),

\[
 \operatorname{tr}M_{\rm an}>
 (b_0-\epsilon_*)p>.005109p,                        \tag{1.13}
\]

and

\[
 {\operatorname{tr}M_{\rm an}\over\|M_{\rm an}\|_{op}}
 \ge {1-\epsilon_*/b_0\over2/500+4\epsilon_*/b_0}
 >19.54.                                            \tag{1.14}
\]

This proves (0.5).

### 1.2 A faithful retained core

Fix

\[
                         \eta_0={1\over2000}.        \tag{1.15}
\]

On \(\{\omega_{\rm an}<\eta_0\}\),

\[
 \eta=\omega_{\rm an}|m|\le\eta_0R.
\]

Since \(\int R\,d\mu=B_0\le p\), deleting this region loses trace at most
\(\eta_0p\).

Let

\[
 M_{\omega}=\int_{\{\omega_{\rm an}\ge\eta_0\}}
                  \eta\theta\theta^T\,d\mu.          \tag{1.16}
\]

Then

\[
 \operatorname{tr}M_{\omega}
 \ge b-\Delta-\eta_0p>.004609p,
\qquad
 {\operatorname{tr}M_{\omega}\over\|M_{\omega}\|_{op}}
 >17.62.                                             \tag{1.17}
\]

On \(\{q<1/2\}\), equation (1.10) gives

\[
 d=R(1-q)>R/2\ge r_G/2,
 \qquad \eta\le r_G<2d.                              \tag{1.18}
\]

Deleting this region loses trace at most \(2\Delta\).  Let

\[
 M_{\mathcal C}=\int_{\mathcal C}\eta\theta\theta^T\,d\mu,
 \qquad
 \mathcal C=\{\omega_{\rm an}\ge\eta_0,\ q\ge1/2\}. \tag{1.19}
\]

Then

\[
 \operatorname{tr}M_{\mathcal C}
 \ge b-3\Delta-\eta_0p
 >(b_0-3\epsilon_*-\eta_0)p
 >.004489p.                                         \tag{1.20}
\]

Since \(M_{\mathcal C}\preceq D\), its operator norm obeys (1.8).  Hence

\[
 {\operatorname{tr}M_{\mathcal C}\over\|M_{\mathcal C}\|_{op}}
 \ge {b_0-3\epsilon_*-\eta_0
       \over (2/500)b_0+4\epsilon_*}
 >17.16.                                            \tag{1.21}
\]

This proves (0.7).

## 2. Pointwise Gaussian-scale regularity

The following elementary bounds require no regularity of the input
function.

**Lemma 2.1 (bounded Gaussian convolution).**  Let
\(\varphi:\mathbb R^n\to\mathbb R^m\) be measurable with
\(|\varphi|\le M\), and put \(a=P_s\varphi\).  Then

\[
 |a|\le M,\qquad
 \|\nabla a\|_{op}\le {M\over\sqrt s},               \tag{2.1}
\]

and for unit vectors \(e,f\) and a unit output covector \(v\),

\[
 \left|D^2_{e,f}\langle v,a\rangle\right|
 \le {\sqrt2M\over s}.                              \tag{2.2}
\]

**Proof.**  Gaussian score differentiation gives

\[
 D_ea(x)={1\over\sqrt s}\mathbb E[Z_e
                 \varphi(x+\sqrt sZ)],              \tag{2.3}
\]

\[
 D^2_{e,f}\langle v,a(x)\rangle
 ={1\over s}\mathbb E[(Z_eZ_f-\langle e,f\rangle)
                  \langle v,\varphi(x+\sqrt sZ)\rangle].          \tag{2.4}
\]

The first bound follows by choosing an output covector attaining the norm
and applying Cauchy--Schwarz.  The centered quadratic Gaussian in (2.4)
has second moment \(1+\langle e,f\rangle^2\le2\).  QED.

Applying Lemma 2.1 with \(M=M_s=I_0/\sqrt s\) gives

\[
 0\le r_G,R\le {I_0\over\sqrt s},\qquad
 |m|\le {I_0\over\sqrt s},                          \tag{2.5}
\]

\[
 |\nabla r_G|,\ |\nabla R|,\ \|\nabla m\|_{op}
 \le {I_0\over s}.                                  \tag{2.6}
\]

On \(\{m\ne0\}\),

\[
 \nabla\theta={P_{\theta^\perp}\nabla m\over|m|},
\qquad
 \|\nabla(\theta\theta^T)\|_{HS}
 =\sqrt2\,\|\nabla\theta\|_{HS}.                    \tag{2.7}
\]

The singular denominator in (2.7) is genuine.  It is removed either by
restricting to \(q\ge1/2\), as in Section 5, or by differentiating the full
matrix density

\[
 Z_{\rm an}(x)=\eta(x)\theta(x)\theta(x)^T
 ={r_G(x)\over R(x)}\,{m(x)m(x)^T\over|m(x)|}.       \tag{2.8}
\]

The map \(m\mapsto mm^T/|m|\), extended by zero, is globally
three-Lipschitz in Hilbert--Schmidt norm.  Thus \(Z_{\rm an}\) is locally
Lipschitz.  Its pointwise derivative also contains
\(\nabla(r_G/R)\), whose ratios are best controlled in Fisher rather than
uniform norm; this is the purpose of Sections 3--4.

## 3. A dimension-free tilted-channel lemma

This section proves (0.9)--(0.10).

Let \(0\le h\le M\), let

\[
 r(x)=P_sh(x)>0,\qquad
 dP_x(z)={h(x+\sqrt sz)\over r(x)}\,d\gamma_n(z),    \tag{3.1}
\]

and let \(u:\mathbb R^n\to\mathbb R^m\) satisfy \(|u|\le1\).  Define

\[
                         v(x)=\mathbb E_{P_x}
                                   u(x+\sqrt sZ).    \tag{3.2}
\]

Put \(D_x=D(P_x\|\gamma_n)\).

**Lemma 3.1 (direction and scalar Fisher bounds).**

\[
 \boxed{
 \int r\|\nabla v\|_{HS}^2d\mu
 \le {8\over s}\int r(1+D_x)d\mu,}                  \tag{3.3}
\]

\[
 \boxed{
 \int {|\nabla r|^2\over r}\,d\mu
 \le {2\over s}\int rD_x\,d\mu.}                    \tag{3.4}
\]

Furthermore, if \(b=\int r\,d\mu\), then

\[
 \boxed{\int rD_x\,d\mu\le b\log{M\over b}.}         \tag{3.5}
\]

**Proof.**  Score differentiation of the quotient (3.2) gives the exact
matrix identity

\[
                         \nabla v
 ={1\over\sqrt s}\operatorname{Cov}_{P_x}(u,Z).      \tag{3.6}
\]

Let \(C=\operatorname{Cov}_{P_x}(u,Z)\).  Matrix Cauchy--Schwarz gives

\[
 \|C\|_{HS}^2
 \le\lambda_{\max}(\operatorname{Cov}_{P_x}Z)\,
       \operatorname{tr}(\operatorname{Cov}_{P_x}u)
 \le\lambda_{\max}(\operatorname{Cov}_{P_x}Z).      \tag{3.7}
\]

For a unit vector \(e\), let
\(\sigma_e^2=\operatorname{Var}_{P_x}\langle e,Z\rangle\).
Data processing and the maximum-entropy property of the Gaussian give

\[
 D_x\ge {1\over2}
 \left(\sigma_e^2-1-\log\sigma_e^2\right).           \tag{3.8}
\]

If \(\sigma_e^2\le4\), it is at most \(4(1+D_x)\).  If it is larger than
four, then
\(\sigma_e^2-1-\log\sigma_e^2\ge\sigma_e^2/4\), so
\(\sigma_e^2\le8D_x\).  Therefore

\[
 \lambda_{\max}(\operatorname{Cov}_{P_x}Z)
 \le8(1+D_x).                                       \tag{3.9}
\]

Equations (3.6)--(3.9) prove (3.3).

Also

\[
                         \nabla\log r
 ={1\over\sqrt s}\mathbb E_{P_x}Z.                  \tag{3.10}
\]

The entropy bound \(D_x\ge|\mathbb E_{P_x}Z|^2/2\) proves (3.4).

Finally,

\[
 \int rD_x\,d\mu
 =\int\operatorname{Ent}_{\gamma_n}
       (h(x+\sqrt s\,\cdot))\,d\mu(x).              \tag{3.11}
\]

The first entropy term is at most \(b\log M\).  Convexity of
\(t\log t\) gives
\(\int r\log r\,d\mu\ge b\log b\).  Subtraction proves (3.5).  QED.

Apply Lemma 3.1 first with

\[
 h=|W|,\qquad u={W\over|W|},\qquad
 r=R,\qquad v={m\over R}.                            \tag{3.12}
\]

This gives (0.9).  Apply it again with

\[
 h=1_G|W|,\qquad u={W\over|W|},\qquad
 r=r_G,\qquad v=v_G={P_s(1_GW)\over r_G}.            \tag{3.13}
\]

This gives (0.10).

## 4. Binary Fisher information of the analytic selector

Fix \(x\), and under the \(|W|\)-tilted channel \(P_x\) from (3.12), let

\[
                         L=1_G(x+\sqrt sZ).
\]

Then

\[
                         P_x(L=1)=\omega_{\rm an}(x).             \tag{4.1}
\]

Let \(P_{x,1},P_{x,0}\) be the two conditional laws of \(Z\), and put
\(\bar z_i=\mathbb E_{P_{x,i}}Z\).  Quotient differentiation gives

\[
 \boxed{
 \nabla\log{\omega_{\rm an}\over1-\omega_{\rm an}}
 ={1\over\sqrt s}(\bar z_1-\bar z_0).}              \tag{4.2}
\]

The variance of the two component means is bounded by their Gaussian
entropy costs:

\[
 \omega(1-\omega)|\bar z_1-\bar z_0|^2
 \le\omega|\bar z_1|^2+(1-\omega)|\bar z_0|^2
\]

\[
 \le2\{\omega D(P_{x,1}\|\gamma_n)
          +(1-\omega)D(P_{x,0}\|\gamma_n)\}.         \tag{4.3}
\]

The mixture entropy chain rule is

\[
 \omega D(P_{x,1}\|\gamma_n)
 +(1-\omega)D(P_{x,0}\|\gamma_n)
 =D(P_x\|\gamma_n)+H_2(\omega),                     \tag{4.4}
\]

where

\[
 H_2(t)=-t\log t-(1-t)\log(1-t)\le\log2.             \tag{4.5}
\]

Multiply (4.2)--(4.4) by \(R\omega(1-\omega)\), integrate, and use
(3.5) with \(h=|W|\).  This proves

\[
 \int R\,{|\nabla\omega_{\rm an}|^2
                 \over\omega_{\rm an}(1-\omega_{\rm an})}\,d\mu
 \le {2B_0\over s}
       \left(\log{M_s\over B_0}+\log2\right).        \tag{4.6}
\]

In particular, the selector cannot change repeatedly without paying binary
Fisher energy.  The estimate is nevertheless an upper budget, not a lower
charge forced by between-level variance.

## 5. Phase energy on the faithful core

Write

\[
 v_0={m\over R}.
\]

On the core \(\mathcal C\), \(|v_0|=q\ge1/2\), and
\(\theta=v_0/|v_0|\).  Hence

\[
 \|\nabla\theta\|_{HS}\le2\|\nabla v_0\|_{HS},\qquad
 \|\nabla(\theta\theta^T)\|_{HS}^2
 \le8\|\nabla v_0\|_{HS}^2.                         \tag{5.1}
\]

Since \(\eta\le R\), (0.9) gives

\[
 \boxed{
 \int_{\mathcal C}\eta
       \|\nabla(\theta\theta^T)\|_{HS}^2d\mu
 \le {64B_0\over s}
       \left(1+\log{M_s\over B_0}\right).}           \tag{5.2}
\]

Using \(B_0\le p\) and (0.12), the right side is at most

\[
                         {C(1+L_\alpha)p\over s}.    \tag{5.3}
\]

The selector transition energy (4.6) has the same scale.

These estimates have the correct Gaussian derivative \(s^{-1}\), contain
no dimension, and contain only the fixed logarithm \(L_\alpha\).  They do
not contain the adverse powers \(\alpha^{-3/2}\) from the earlier scalar
angular-stability audit.

## 6. Between-level variance transfers to a smooth heat direction

Let

\[
 a=\int\eta\,d\mu,\qquad d\nu={\eta\over a}\,d\mu,\qquad
 A(x)=\theta(x)\theta(x)^T,\qquad C(x)=v_G(x)v_G(x)^T.             \tag{6.1}
\]

The conditional good-direction mean satisfies the pointwise alignment
identity

\[
 r_G(1-\langle\theta,v_G\rangle)
 =P_s\!\left[1_G|W|(1-\langle\theta,u\rangle)\right]
 \le d.                                                \tag{6.2}
\]

Since \(|v_G|\le1\),

\[
 r_G|v_G-\theta|^2\le2d.                             \tag{6.3}
\]

Also

\[
 \|A-C\|_{HS}
 \le(|\theta|+|v_G|)|\theta-v_G|
 \le2|\theta-v_G|.                                  \tag{6.4}
\]

Because \(\eta\le r_G\), integration gives

\[
                         \mathbb E_\nu\|A-C\|_{HS}^2
 \le {8\Delta\over a}.                              \tag{6.5}
\]

Let \(\mathsf P\) denote conditional expectation onto
\(\sigma(F_0)\), followed by removal of the mean.  It is an orthogonal
projection in the matrix-valued \(L^2(\nu)\) space.  Therefore

\[
 \left|
 \|\mathsf PA\|_{L^2(\nu)}-\|\mathsf PC\|_{L^2(\nu)}
 \right|
 \le\sqrt{8\Delta/a}.                               \tag{6.6}
\]

Suppose the between-level variance of the analytic physical normal matrix
obeys

\[
                         \|\mathsf PA\|_{L^2(\nu)}^2
 \ge {8\over17}.                                    \tag{6.7}
\]

Using \(a\ge(b_0-\epsilon_*)p\) and
\(\Delta\le\epsilon_*p\), equations (6.6)--(6.7) give

\[
 \|\mathsf PC\|_{L^2(\nu)}^2
 \ge\left(
 \sqrt{8/17}
 -\sqrt{8\epsilon_*/(b_0-\epsilon_*)}
 \right)^2
 >0.143.                                            \tag{6.8}
\]

This proves (0.13).

### 6.1 The unavoidable capacity loss

Equation (6.8) is a fixed variance statement, while (5.2) is an upper
Dirichlet budget.  A lower Dirichlet charge would require an inequality of
the form

\[
 a\,\operatorname{Var}_\nu H
 \le C_{\rm cap}\int_{\mathcal C}\eta|\nabla H|^2d\mu             \tag{6.9}
\]

for the selected matrix field.  The measure \(\nu\) is not log-concave;
it is a boundary-flux tilt with two additional threshold restrictions.
Neither Gaussian convolution nor (4.6) proves (6.9).

Even granting the optimistic ambient-scale value
\(C_{\rm cap}=CK\), (6.8) yields only

\[
 \int_{\mathcal C}\eta|\nabla H|^2d\mu
 \ge {c\,a\over K}\ge {c'p\over K}
 ={c'\alpha p\over s}.                              \tag{6.10}
\]

The natural scale-free Gaussian charge is \(s\) times this energy, namely

\[
                         s\int\eta|\nabla H|^2
 \ge c'\alpha p=10^{-10}c'p.                        \tag{6.11}
\]

This is six decimal orders below \(10^{-4}p\).  An argument which invokes
the ambient Poincare inequality here has therefore reintroduced the old
power obstruction.  Proving (6.9) with \(C_{\rm cap}=O(s)\) would avoid the
loss, but such a selected-measure capacity statement is a new
extremality theorem, not a consequence of the Gaussian derivative bounds.

There is an additional level-space obstruction.  The coarea matrices are
the disintegration of \(\eta\theta\theta^T\mu\) under \(F_0\).  Differentiating
them in the level variable requires the vector field

\[
 {\eta\over|m|}\,(\theta^TH\theta)\theta
 =\omega_{\rm an}(\theta^TH\theta)\theta,            \tag{6.12}
\]

for matrix test directions \(H\).  Although its magnitude is bounded, its
weighted divergence contains the full weighted mean curvature, focal
terms, and critical values of \(F_0\).  Bounds (0.9)--(0.11) control none
of those divergences.  Thus a level-BV conversion is exactly another form
of the missing physical incidence theorem.

## 7. An exact heat-generated local phase profile

This section shows that the information estimates alone cannot force a
Bernstein charge.  The model is physical and heat-generated, not an
abstract graph, but it is deliberately not asserted to be a high-\(K\)
near-Cheeger counterexample.

### 7.1 One log-affine planar facet

Work in a product chart with signed normal coordinate \(t\), and suppose

\[
 S=\{t\le0\},\qquad
 e^{-V(z,t)}=a(z)e^{-\kappa t}.                     \tag{7.1}
\]

Ignore the remote boundary of the chart; all formulas below are exact for
the full log-affine model and have error \(O(e^{-cL^2})\) in a chart whose
facet core has Gaussian reach \(L\sqrt s\).

Given an observation normal coordinate \(y\), completing the square shows
that the posterior normal coordinate is

\[
                         t\mid y\sim N(y-s\kappa,s).               \tag{7.2}
\]

Therefore

\[
 g(y)=\Phi\!\left({s\kappa-y\over\sqrt s}\right),\qquad
 z(y)={s\kappa-y\over\sqrt s},                     \tag{7.3}
\]

\[
 W(y)=-{I(g(y))\over\sqrt s}\,n,\qquad
 e=s|\nabla z|^2=1,\qquad
 \nabla^2z=0.                                       \tag{7.4}
\]

The physical smoothed label is

\[
 F_0(x)=\mathbb E g(x+\sqrt sZ)
 =\Phi\!\left({s\kappa-x\over\sqrt{2s}}\right),      \tag{7.5}
\]

and

\[
 m(x)=-{1\over\sqrt{2s}}\,
 \varphi\!\left({s\kappa-x\over\sqrt{2s}}\right)n,
\qquad R(x)=|m(x)|,\qquad d(x)=0.                   \tag{7.6}
\]

For the central good set \(G\), the restricted flux \(r_G\) is a universal
function of the profile coordinate in (7.5), independent of \(\kappa\).
The analytic selector is

\[
                         \omega_{\rm an}={r_G\over R},             \tag{7.7}
\]

and the selected normal is the constant line \(\mathbb Rn\).

On the physical level \(F_0=r\), write \(z_r=\Phi^{-1}(r)\).  Equation
(7.5) gives

\[
                         x=s\kappa-\sqrt{2s}\,z_r.   \tag{7.8}
\]

The ambient density on this level contributes the factor

\[
 e^{-\kappa x}
 =\exp\{-s\kappa^2+\sqrt{2s}\kappa z_r\}.            \tag{7.9}
\]

### 7.2 Many phases

Take disjoint planar facet cores with unit normals \(n_i\), log-affine
normal slopes \(\kappa_i\), and weighted tangential coefficients \(A_i\).
Up to the common universal profile factor supplied by \(r_G\), their
physical level matrix is

\[
 M_r=\sum_i
 A_i\exp\{-s\kappa_i^2+\sqrt{2s}\kappa_i z_r\}
 n_in_i^T.                                         \tag{7.10}
\]

Put

\[
                         c_i=\sqrt s\,\kappa_i.
\]

Then the normalized phase weights form the softmax

\[
 \pi_i(z)=
 {A_i e^{-c_i^2+\sqrt2c_i z}
  \over\sum_jA_j e^{-c_j^2+\sqrt2c_j z}}.           \tag{7.11}
\]

Choose eighteen well-conditioned normals, choose separated slopes \(c_i\),
and tune the positive coefficients \(A_i\) so that consecutive phases cross
at prescribed profile levels.  As the slope separations tend to infinity,
one phase dominates away from arbitrarily narrow crossing intervals.  Hence

\[
 \mathsf W_{\rm same}\longrightarrow0,\qquad
 \mathsf B_{\rm between}\longrightarrow
 1-\operatorname{tr}\!\left(
 {1\over18}\sum_{i=1}^{18}n_in_i^T\right)^2,         \tag{7.12}
\]

which is \(17/18\) for orthonormal normals.

Every facet core in this construction has the exact heat properties

\[
                         e=1,\qquad \nabla^2z=0,\qquad d=0.        \tag{7.13}
\]

Thus neither the Hessian square nor the eikonal square in the Bernstein
identity charges the cross-level rotation on the cores.  For a genuinely
log-affine potential, the curvature square also vanishes there.

The model can be embedded, with \(O(e^{-cL^2})\) error, in a log-concave
probability \(e^{-\langle b,x\rangle}1_Cdx\) by taking a large convex
support \(C\), choosing facet cores at distance \(L\sqrt s\) from its
boundary and from their ridges, and arranging
\(\kappa_i=\langle b,n_i\rangle\).  The set \(S\) is completed away from
the cores by finite-perimeter patches.  Gaussian observations which stay
inside a core see exactly (7.2)--(7.7).

The completion necessarily creates ridges, ends, or transition regions.
Those are the only locations where a perimeter saving or a singular
turning charge can occur.  This is precisely the physical datum absent from
the Gaussian derivative and information inequalities.

The example is not a near-Cheeger set and has no claim of arbitrarily large
isotropic Poincare constant.  Its conclusion is narrower and exact:
heat generation, analyticity of \(\omega_{\rm an}\), centrality, and
near-eikonal posterior states do not by themselves turn between-level
projector variance into a Bernstein charge.  Global near-extremality must
force and then exploit the ridge or transition charge.

## 8. Exact remaining statement

The selector regularity problem is now reduced to the following
extremality-dependent alternative.

1. If the physical rank is paid by same-level dispersion, prove a
   bounded-overlap ridge, focal, or short-bridge conductance estimate and
   apply the finite bevel lemma.
2. If it is paid by between-level dispersion, use (0.13) together with
   near-Cheeger extremality to prove a selected-measure capacity bound with
   scale \(O(s)\), or prove a full level-BV/focal charge.  The ambient scale
   \(K=s/\alpha\) is insufficient by (6.11).

The analytic selector eliminates arbitrary hard-threshold rotation and
retains more rank than the original construction.  The binary Fisher
identity quantifies its spatial motion with no dimension loss.  What
remains is not regularity of the selector; it is the physical capacity of
the phase transition.  The log-affine model shows that this capacity lives
at the incidence skeleton, exactly where the existing finite bevel
construction operates.

For the blocked-list audit, the conservative numerical certificate used
throughout is

\[
 \boxed{
 \begin{array}{c|c|c}
 \text{physical selector} & \text{trace lower bound}
                           & \text{effective-rank lower bound}\\ \hline
 \omega_{\rm an}=r_G/R
       & .005109\,p & 19.54\\
 \omega_{\rm an}\,1_{\{\omega_{\rm an}\ge1/2000\}}
       & .004609\,p & 17.62\\
 \omega_{\rm an}\,1_{\{\omega_{\rm an}\ge1/2000,\ |m|/R\ge1/2\}}
       & .004489\,p & 17.16
 \end{array}}                                                   \tag{8.1}
\]

These decimals use only \(b_0>.00517\) and
\(\Delta<6.02\cdot10^{-5}p\), not the slightly stronger decimal obtained
by evaluating the exact formula for \(b_0\).
