# Exact tests of the final-time alignment inequality

## 1. Setup and verdict

For the standard stochastic-localization family

\[
 d\mu_{t,\theta}(x)=Z(t,\theta)^{-1}
 \exp\{\langle\theta,x\rangle-t|x|^2/2\}\,d\mu(x),
\]

let \(\theta_t\) solve

\[
 d\theta_t=dB_t+a(t,\theta_t)\,dt,
 \qquad a(t,\theta)=\mathbb E_{t,\theta}X,
 \qquad \theta_0=0.
\]

Write \(A_t=\operatorname{Cov}_{\mu_{t,\theta_t}}(X)\), and let

\[
 M_0=I,\qquad \dot M_t=A_tM_t.
\]

For a fixed function \(f\), put

\[
 c_t=\operatorname{Cov}_{\mu_{t,\theta_t}}(X,f).
\]

The proposed final-time alignment estimate is

\[
 \mathbb E\left[
 \frac{\langle M_1b,c_1\rangle^2}{|c_1|^2}
 \mathbf 1_{\{c_1\ne0\}}
 \right]\le C. \tag{FA}
\]

No counterexample was found.  The calculations below prove (FA), without
using an ambient trace estimate, in three nontrivial situations:

1. every Gaussian model (with the exact dependence on covariance);
2. an isotropic interval times an arbitrary number of Gaussian factors,
   for its genuine first eigenfunction, which is almost affine;
3. arbitrary tensor products of uniformly bounded-diameter factors,
   including products of fixed-dimensional regular simplexes and triangular
   wedges, for every \(f\) and every \(b\).

They also identify an important limitation.  A fixed-time version of (FA)
is not affine invariant: an anisotropic Gaussian has alignment quotient
\((1+\sigma)^2\) in a covariance eigen-direction of variance \(\sigma\).
Thus isotropic normalization (or a covariance-scaled localization time) is
essential.

The calculations do **not** prove (FA) for one regular simplex whose
dimension tends to infinity.  The usual symmetry-plus-trace verification
of that example is valid, but it is not an independent directional
mechanism.

## 2. Two deterministic comparison lemmas

### Lemma 2.1 (posterior covariance ceiling)

Suppose that, along a localization path,

\[
 0\preceq A_t\preceq K(t)I\qquad(0\le t\le1)
\]

with \(K\) integrable.  Then, pathwise, for every unit vector \(b\),

\[
 |M_1b|^2\le
 \exp\left(2\int_0^1K(t)\,dt\right). \tag{2.1}
\]

Consequently, for every \(f\),

\[
 \frac{\langle M_1b,c_1\rangle^2}{|c_1|^2}
 \mathbf 1_{\{c_1\ne0\}}
 \le \exp\left(2\int_0^1K(t)\,dt\right). \tag{2.2}
\]

**Proof.**  Set \(w_t=M_tb\).  Since \(A_t\) is symmetric and positive
semidefinite,

\[
 \frac d{dt}|w_t|^2=2\langle w_t,A_tw_t\rangle
 \le2K(t)|w_t|^2.
\]

Gronwall gives (2.1).  Cauchy--Schwarz in Euclidean space gives

\[
 \langle M_1b,c_1\rangle^2/|c_1|^2\le |M_1b|^2
\]

when \(c_1\ne0\), proving (2.2).  No commutativity of the matrices
\(A_t\) is used. \(\square\)

### Lemma 2.2 (bounded support)

If the support of a probability measure has Euclidean diameter at most
\(D\), then every exponential-quadratic tilt of it has covariance at most

\[
 \operatorname{Cov}(X)\preceq \frac{D^2}{4}I. \tag{2.3}
\]

Hence its standard localization flow satisfies, pathwise,

\[
 |M_1b|^2\le e^{D^2/2}. \tag{2.4}
\]

**Proof.**  For every unit \(u\), the range of \(\langle u,X\rangle\)
has length at most \(D\).  Popoviciu's inequality gives
\(\operatorname{Var}\langle u,X\rangle\le D^2/4\).  This remains true
after every tilt because the support is unchanged.  Apply Lemma 2.1 with
\(K=D^2/4\). \(\square\)

A further immediate class, useful as a check on the mechanism, consists of
\(\rho\)-strongly log-concave measures.  If \(\nabla^2V\succeq\rho I\),
then the posterior potential has Hessian at least \((\rho+t)I\), so
Brascamp--Lieb gives \(A_t\preceq(\rho+t)^{-1}I\).  Lemma 2.1 yields

\[
 |M_1b|^2\le\left(\frac{\rho+1}{\rho}\right)^2. \tag{2.5}
\]

In particular, (FA) holds dimension-freely for the non-Gaussian class
\(\rho\ge\rho_0>0\), without any eigenfunction hypothesis.

## 3. Exact Gaussian calculation

Let \(\mu=N(0,\Sigma)\), where \(\Sigma\) is positive definite.  Completing
the square gives

\[
 A_t=(\Sigma^{-1}+tI)^{-1}.
\]

This is deterministic.  Moreover

\[
 M_t=I+t\Sigma, \tag{3.1}
\]

because

\[
 (\Sigma^{-1}+tI)^{-1}(I+t\Sigma)=\Sigma.
\]

Let \(v\) be a unit eigenvector of \(\Sigma\) with eigenvalue \(\sigma\),
and normalize the linear eigenfunction by

\[
 f(x)=\frac{\langle v,x\rangle}{\sqrt\sigma}.
\]

Then

\[
 c_1=A_1\frac v{\sqrt\sigma}
 =\frac{\sqrt\sigma}{1+\sigma}v,
 \qquad M_1v=(1+\sigma)v.
\]

Therefore the alignment quotient is exactly

\[
 \frac{\langle M_1v,c_1\rangle^2}{|c_1|^2}
 =(1+\sigma)^2. \tag{3.2}
\]

For an isotropic Gaussian, \(\sigma=1\), so the value is exactly \(4\).
Equation (3.2) is also a sharp warning: at localization time one the same
claim cannot hold uniformly before isotropic normalization.

For the isotropic Gaussian, in fact \(M_1=2I\), and hence the quotient is
at most \(4\) for **every** function \(f\), not just for a linear
eigenfunction.

## 4. A genuine near-linear first eigenfunction

Let \(L=\sqrt3\), let \(\nu\) be uniform on \([-L,L]\), and set

\[
 \mu=\nu\otimes\gamma_{n-1}.
\]

Both factors are isotropic, so \(\mu\) is isotropic.  Put

\[
 k=\frac\pi{2L},\qquad
 f(x)=\sqrt2\sin(kx_1). \tag{4.1}
\]

The Neumann spectrum of \([-L,L]\) starts at \(k^2=\pi^2/12\), while the
Gaussian gap is \(1\).  Tensorization of the generators therefore shows
that (4.1) is a normalized genuine first eigenfunction of \(\mu\), with

\[
 \lambda=\frac{\pi^2}{12}. \tag{4.2}
\]

Its linear projection is unusually large.  Since
\(\mathbb E_\nu X_1^2=1\), direct integration gives

\[
 a:=\mathbb E[X_1f(X)]
 =\frac{\sqrt2}{L}\int_0^Lx\sin(kx)\,dx
 =\frac{4\sqrt6}{\pi^2}. \tag{4.3}
\]

Thus, for the orthogonal decomposition \(f=ax_1+r\),

\[
 \|r\|_2^2=1-a^2=1-\frac{96}{\pi^4}<0.015. \tag{4.4}
\]

This is a genuine low-mode test in precisely the near-linear branch for
which (FA) was proposed.

Localization preserves the product decomposition.  Write \(v_t\) for the
posterior variance of the interval coordinate and

\[
 m_t=\exp\left(\int_0^t v_s\,ds\right).
\]

Then

\[
 A_t=\operatorname{diag}\left(v_t,(1+t)^{-1}I_{n-1}\right),
 \qquad
 M_t=\operatorname{diag}\left(m_t,(1+t)I_{n-1}\right). \tag{4.5}
\]

Every posterior of the first coordinate is supported on an interval of
length \(2\sqrt3\), so Popoviciu gives \(v_t\le3\) and therefore

\[
 1\le m_1\le e^3. \tag{4.6}
\]

Because \(f\) depends only on \(x_1\), posterior product factorization
gives

\[
 c_t=(\operatorname{Cov}_{\nu_t}(X_1,f),0,\ldots,0). \tag{4.7}
\]

The first component is strictly positive: for independent
\(U,U'\sim\nu_t\),

\[
 2\operatorname{Cov}(U,f(U))
 =\mathbb E[(U-U')(f(U)-f(U'))]>0, \tag{4.8}
\]

because the posterior is nondegenerate and \(f\) is strictly increasing
on the interior of \([-L,L]\).  Taking \(b=e_1\), (4.5)--(4.8) give the
exact pathwise identity

\[
 \frac{\langle M_1b,c_1\rangle^2}{|c_1|^2}=m_1^2\le e^6. \tag{4.9}
\]

This proves (FA) with a universal constant for this genuine first mode in
every dimension.  The calculation uses neither a trace estimate nor
averaging over rotations.

## 5. Product simplexes and triangular wedges

Let \(\mu=\mu_1\otimes\cdots\otimes\mu_N\) on an orthogonal decomposition
\(E_1\oplus\cdots\oplus E_N\).  Standard localization factorizes:

\[
 A_t=\operatorname{diag}(A_{1,t},\ldots,A_{N,t}),
 \qquad
 M_t=\operatorname{diag}(M_{1,t},\ldots,M_{N,t}). \tag{5.1}
\]

If every factor has support diameter at most \(D\), Lemma 2.2 applied
blockwise gives

\[
 \|M_1\|_{\mathrm{op}}^2
 =\max_j\|M_{j,1}\|_{\mathrm{op}}^2
 \le e^{D^2/2}. \tag{5.2}
\]

Consequently (FA) holds pathwise, with the same bound, for every function
\(f\), every deterministic \(b\), and every number of factors.  This is a
dimension-free high-dimensional class although its irreducible blocks have
bounded dimension.

For a concrete simplex calculation, let \(K_d\) be a regular
\(d\)-simplex in isotropic position and let \(v_0,\ldots,v_d\) be its
vertices.  The Dirichlet second-moment formula for uniform barycentric
coordinates gives

\[
 \operatorname{Cov}(X)
 =\frac1{(d+1)(d+2)}\sum_{i=0}^d v_iv_i^T. \tag{5.3}
\]

For a centered regular simplex, isotropy in (5.3) implies

\[
 |v_i|^2=d(d+2),\qquad
 |v_i-v_j|^2=2(d+1)(d+2)\quad(i\ne j). \tag{5.4}
\]

Thus

\[
 \operatorname{diam}(K_d)^2=2(d+1)(d+2). \tag{5.5}
\]

For an arbitrary tensor product of copies of this fixed simplex,
(5.2) becomes

\[
 \frac{\langle M_1b,c_1\rangle^2}{|c_1|^2}
 \mathbf1_{\{c_1\ne0\}}
 \le \exp\{(d+1)(d+2)\}. \tag{5.6}
\]

The constant is independent of the number of factors.  When \(d=2\), an
isotropic triangular wedge has side length \(2\sqrt6\), and (5.6) reads

\[
 \frac{\langle M_1b,c_1\rangle^2}{|c_1|^2}
 \mathbf1_{\{c_1\ne0\}}\le e^{12}. \tag{5.7}
\]

Again this holds for every \(f\), including all genuine or approximate
first modes of the tensor product.

The same argument for one simplex with \(d=n\) gives only
\(e^{(n+1)(n+2)}\), so it supplies no dimension-free information there.
The dimension-free regular-simplex check in the main research note instead
uses irreducible symmetry to turn an ambient trace estimate into a
directional estimate.  That check is consistent, but it should not be
mistaken for a direct control of the alignment process.

## 6. Consequences for the search for a counterexample

These tests impose four concrete requirements on any counterexample to
(FA).

1. It cannot be generated by a uniformly strongly log-concave posterior
   ceiling; somewhere along typical amplified paths the top posterior
   covariance must exceed every fixed constant.
2. It cannot be a tensor product of uniformly bounded-diameter blocks,
   even if the number of blocks tends to infinity.
3. Exact Gaussian directions and the almost-affine interval first mode do
   not create the required misalignment.  A counterexample must exploit
   genuinely rotating covariance eigenspaces or a single irreducible block
   whose dimension grows.
4. High-dimensional regular simplex symmetry conceals the selected
   direction through averaging.  An informative wedge/cone candidate must
   break that irreducible symmetry while retaining one macroscopic
   near-linear first mode.

The natural next explicit model is therefore an axially symmetric,
isotropic log-concave cone with growing irreducible dimension, rather than
a product wedge.  One tractable candidate has unnormalized density

\[
 \mathbf1_{\{s>0\}}
 \exp\left\{-s-\frac{|z|^2}{2s}\right\}\,ds\,dz,
 \qquad z\in\mathbb R^{n-1}. \tag{6.1}
\]

The potential in (6.1) is convex (the quadratic-over-linear term is a
perspective), and conditionally

\[
 Z\mid S=s\sim N(0,sI),\qquad
 S\sim\operatorname{Gamma}\left(\frac{n+1}{2},1\right). \tag{6.2}
\]

Put

\[
 \alpha=\frac{n+1}{2},\qquad
 X_0=\frac{S-\alpha}{\sqrt\alpha},\qquad
 X_z=\frac Z{\sqrt\alpha}.
\]

Then \(X=(X_0,X_z)\) is exactly isotropic.  The localized posterior can be
reduced to one-dimensional quadrature.  At time \(t\) and natural
parameter \(\theta=(u,v)\in\mathbb R\times\mathbb R^{n-1}\), define

\[
 q_t(s)=\frac{\alpha s}{\alpha+ts}. \tag{6.3}
\]

Conditionally on \(S=s\), the posterior law is

\[
 Z\mid(S=s,t,u,v)
 \sim N\left(\frac{q_t(s)}{\sqrt\alpha}v,
                   q_t(s)I_{n-1}\right). \tag{6.4}
\]

The one-dimensional posterior density of \(S\), up to normalization, is

\[
 \begin{split}
 \pi_{t,u,v}(s)\propto{}&s^{\alpha-1}e^{-s}
 \left(1+\frac{ts}{\alpha}\right)^{-(\alpha-1)}\\
 &\times\exp\left\{
 \frac{u(s-\alpha)}{\sqrt\alpha}
 -\frac{t(s-\alpha)^2}{2\alpha}
 +\frac{q_t(s)|v|^2}{2\alpha}
 \right\}\mathbf1_{\{s>0\}}. \tag{6.5}
 \end{split}
\]

All entries of the drift and covariance now follow from scalar moments
under (6.5).  If \(\mathbb E_*\), \(\operatorname{Var}_*\), and
\(\operatorname{Cov}_*\) denote moments under (6.5), and

\[
 h=\frac{\mathbb E_*q_t(S)}{\alpha},
\]

then

\[
 a(t,u,v)=
 \begin{pmatrix}
 (\mathbb E_*S-\alpha)/\sqrt\alpha\\ hv
 \end{pmatrix}, \tag{6.6}
\]

and the covariance has the exact arrowhead form

\[
 A(t,u,v)=
 \begin{pmatrix}
 \operatorname{Var}_*(S)/\alpha &
 \dfrac{\operatorname{Cov}_*(S,q_t(S))}{\alpha\sqrt\alpha}v^T\\[6pt]
 \dfrac{\operatorname{Cov}_*(S,q_t(S))}{\alpha\sqrt\alpha}v &
 hI+\dfrac{\operatorname{Var}_*(q_t(S))}{\alpha^2}vv^T
 \end{pmatrix}. \tag{6.7}
\]

Thus the exact parameter process is

\[
 du_t=dB_t^{(0)}+\frac{\mathbb E_*S-\alpha}{\sqrt\alpha}\,dt,
 \qquad
 dv_t=dB_t^{(z)}+h_tv_t\,dt, \tag{6.8}
\]

coupled to \(\dot M_t=A(t,u_t,v_t)M_t\).

There is also an exact curvature check on the scalar reduction.  If
\(U=-\log\pi_{t,u,v}\) (up to an additive constant), differentiation of
(6.5) gives

\[
 U''(s)=\frac t\alpha
 +(\alpha-1)\left\{\frac1{s^2}
 -\frac{t^2}{(\alpha+ts)^2}\right\}
 +\frac{\alpha t|v|^2}{(\alpha+ts)^3}
 \ge\frac t\alpha. \tag{6.9}
\]

In particular \(\operatorname{Var}_*(S)\le\alpha/t\), reproducing the
standard axial covariance ceiling but not a constant one.  Formula (6.7)
also displays the actual difficulty: the axial direction couples to the
instantaneous random direction \(v_t/|v_t|\), which rotates along the
path.  Consequently the matrices \(A_t\) do not commute and the product
integral \(M_t\) is not an endpoint function.  Establishing or refuting
(FA) for (6.1) requires control of this axial/radial path coupling;
fixed-endpoint posterior formulas alone do not determine the required
conditional second moment.

The rotating-column problem nevertheless has an exact scalar Markov
reduction.  Abbreviate the four covariance coefficients in (6.7) by

\[
 \sigma=\frac{\operatorname{Var}_*(S)}\alpha,
 \qquad
 \eta=\frac{\operatorname{Cov}_*(S,q_t(S))}{\alpha\sqrt\alpha},
 \qquad
 h=\frac{\mathbb E_*q_t(S)}\alpha,
 \qquad
 k=\frac{\operatorname{Var}_*(q_t(S))}{\alpha^2}. \tag{6.10}
\]

Let the axial column of the product integral be

\[
 M_te_0=(p_t,z_t)\in\mathbb R\times\mathbb R^{n-1}
\]

and set

\[
 R_t=|v_t|^2,\qquad Q_t=\langle v_t,z_t\rangle,
 \qquad Z_t=|z_t|^2. \tag{6.11}
\]

Writing \(d=n-1=2(\alpha-1)\), direct differentiation of
\(\dot M=AM\) and It\^o differentiation of \(R,Q\) give

\[
 \begin{aligned}
 du={}&dB^{(0)}+\frac{\mathbb E_*S-\alpha}{\sqrt\alpha}\,dt,\\
 dR={}&2\langle v,dB^{(z)}\rangle+(d+2hR)\,dt,\\
 dp={}&(\sigma p+\eta Q)\,dt,\\
 dQ={}&\langle z,dB^{(z)}\rangle
       +(2hQ+\eta Rp+kRQ)\,dt,\\
 dZ={}&(2\eta pQ+2hZ+2kQ^2)\,dt. \tag{6.12}
 \end{aligned}
\]

The initial data are \((u,R,p,Q,Z)=(0,0,1,0,0)\).  The two transverse
noise terms in (6.12) have brackets

\[
 d\langle R\rangle_t=4R_tdt,
 \qquad d\langle Q\rangle_t=Z_tdt,
 \qquad d\langle R,Q\rangle_t=2Q_tdt. \tag{6.13}
\]

Thus (6.12)--(6.13), together with the independent axial Brownian motion,
are a closed finite-dimensional scalar diffusion; the growing ambient
dimension appears only through \(d=2(\alpha-1)\) and through the scalar
posterior coefficients.

For the axial linear test \(f(X)=X_0\), the posterior gradient is

\[
 c_t=A_te_0=(\sigma_t,\eta_tv_t),
\]

so its alignment quotient has the exact form

\[
 \boxed{
 \frac{\langle M_te_0,c_t\rangle^2}{|c_t|^2}
 =\frac{(\sigma_tp_t+\eta_tQ_t)^2}
        {\sigma_t^2+\eta_t^2R_t}.} \tag{6.14}
\]

In particular the numerator is \((\dot p_t)^2\).  Formula (6.14) turns
the noncommutative cone test into a concrete scalar SDE estimate, without
discarding the rotational noise.  The elementary coefficient bounds

\[
 0\le\eta\le\frac\sigma{\sqrt\alpha},
 \qquad 0\le k\le\frac\sigma\alpha \tag{6.15}
\]

follow because \(q_t(s)\) is increasing and 1-Lipschitz, while Cauchy--
Schwarz gives the corresponding covariance bound.  They show that
rotation is coupled to the natural ratios \(R/\alpha\) and
\(Z/\alpha\), not directly to the ambient dimension.  They do not by
themselves bound (6.14), since the scalar axial variance \(\sigma_t\) can
still have rare spikes of order \(1/t\).

There is a sharper reduction which discards only angular energy that the
alignment quotient cannot see.  For \(R>0\), put

\[
 P=\frac{Q^2}{R},\qquad H=Z-P\ge0. \tag{6.16}
\]

Cauchy--Schwarz in \(\mathbb R^2\), applied to (6.14), gives the pathwise
bound

\[
 \frac{(\sigma p+\eta Q)^2}{\sigma^2+\eta^2R}
 \le p^2+P. \tag{6.17}
\]

Thus the potentially large angular part \(H\) of \(M_te_0\) is irrelevant
except insofar as the Brownian rotation can feed it back into \(P\).
Applying It\^o to \(P=Q^2/R\) in (6.12)--(6.13) gives, up to a local
martingale,

\[
 \begin{aligned}
 \frac d{dt}P={}&2\eta pQ+2hP+2kQ^2
               +\frac{H-(d-1)P}{R},\\
 \frac d{dt}H={}&2hH-\frac{H-(d-1)P}{R},\\
 \frac d{dt}p^2={}&2\sigma p^2+2\eta pQ. \tag{6.18}
 \end{aligned}
\]

For rigor one first stops at \(R=\varepsilon\) and bounded state, and
then lets \(\varepsilon\downarrow0\); when \(d\ge2\), \(R_t>0\) for every
\(t>0\) almost surely.  The sole rotation-exchange term is

\[
 \frac{H-(d-1)P}{R}. \tag{6.19}
\]

It enters the visible radial energy \(P\) and the invisible angular energy
\(H\) with opposite signs.  Consequently a proof for this cone can be
reduced further to a dimension-free Lyapunov bound for \(p^2+P\) which
charges positive excursions of (6.19) to the subsequent damping in the
second equation of (6.18).  This is substantially narrower than bounding
the full column norm \(p^2+Z\), and it does not invoke an ambient trace
estimate.

## 7. Bottom line

The final-time alignment inequality survives exact tests on a genuine
near-linear first eigenfunction and on noncommuting bounded-block product
flows.  These tests are stronger than a trace-only consistency check, but
they do not address the only plausible failure regime: a growing,
irreducible, symmetry-broken block with a covariance spike that rotates
into the direction selected by the low mode.  The cone model (6.1) is a
precise next stress test.
