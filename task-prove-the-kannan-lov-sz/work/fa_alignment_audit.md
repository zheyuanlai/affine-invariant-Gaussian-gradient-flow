# Audit of the final-time localization alignment target

## 1. Framework and a fully justified regular case

Let \(\mu\) be a centered probability measure on \(\mathbb R^n\), and set

\[
 Z(t,\theta)=\int \exp\!\left(\langle\theta,x\rangle-
                    \frac t2|x|^2\right)\,d\mu(x),
 \qquad
 \mu_{t,\theta}(dx)=Z(t,\theta)^{-1}
 e^{\langle\theta,x\rangle-t|x|^2/2}\mu(dx).
 \tag{1.1}
\]

Write

\[
 a(t,\theta)=\mathbb E_{t,\theta}X=\nabla_\theta\log Z(t,\theta),
 \qquad
 A(t,\theta)=\operatorname{Cov}_{t,\theta}(X)
             =D_\theta^2\log Z(t,\theta).
 \tag{1.2}
\]

On a filtered probability space carrying standard \(n\)-dimensional
Brownian motion, standard stochastic localization is

\[
 d\theta_t=dB_t+a(t,\theta_t)\,dt,\qquad \theta_0=0,
 \qquad \mu_t=\mu_{t,\theta_t}.
 \tag{1.3}
\]

The cleanest completely elementary regularity regime is:

* \(\mu\) has compact support contained in a ball of radius \(R\);
* \(f\) is bounded and Borel (smoothness is only needed later when the
  original generator is invoked).

Then \(a\) is globally bounded and globally Lipschitz in \(\theta\), since
\(\|D_\theta a\|_{\mathrm{op}}=\|A\|_{\mathrm{op}}\le R^2\).  Hence
(1.3) has a unique nonexplosive strong solution.  All posterior moments
below are bounded, and every stochastic integral displayed below is a true
square-integrable martingale.  This regime is sufficient to derive all
identities without hidden stopping or integrability assumptions.

For a general log-concave \(\mu\), (1.3) is the usual Eldan process.  For
every \(\varepsilon>0\), the posterior is \(t\)-strongly log-concave at
time \(t>0\), and therefore

\[
 0\preceq A(t,\theta)\preceq t^{-1}I.
 \tag{1.4}
\]

Thus all coefficients are smooth and locally Lipschitz on
\([\varepsilon,1]\times\mathbb R^n\).  The calculations below remain
valid as stopped local-martingale identities if

\[
 \tau_R=\inf\{t:\ |\theta_t|+\|M_t\|_{\mathrm{op}}
                +|c_t|+\|D_t\|_{\mathrm{HS}}\ge R\}\wedge1.
 \tag{1.5}
\]

Passing from \(t\wedge\tau_R\) to \(t\) requires uniform integrability;
it is not a formal consequence of local It\^o calculus.  One rigorous
route is first to prove the desired estimate uniformly for compact
truncations and bounded truncations of \(f\), and then use local uniform
convergence of the partition functions and their derivatives on
\([\varepsilon,1]\times\{\theta:|\theta|\le L\}\), followed by
\(\varepsilon\downarrow0\), \(L\uparrow\infty\).  Any proposed proof of
(FA) in full generality has to include this uniform-integrability step.

## 2. Posterior derivatives and the martingale

For fixed \(f\in L^2(\mu)\), define

\[
 F(t,\theta)=\mathbb E_{t,\theta}f,
 \qquad
 c(t,\theta)=\operatorname{Cov}_{t,\theta}(X,f),
 \tag{2.1}
\]

and along the process write \(F_t=F(t,\theta_t)\), \(c_t=c(t,\theta_t)\).
Exponential-family differentiation gives the two exact identities

\[
 \nabla_\theta F=c,
 \qquad
 D_\theta^2F=D,
 \quad
 D=\mathbb E_{t,\theta}
 \big[(f-F)(X-a)(X-a)^T\big].
 \tag{2.2}
\]

Both the numerator and denominator defining \(F\) solve the backward heat
equation.  Consequently

\[
 \partial_tF+\frac12\Delta_\theta F+a\cdot\nabla_\theta F=0.
 \tag{2.3}
\]

It\^o's formula, followed by differentiation of (2.3), yields

\[
 dF_t=c_t\cdot dB_t,
 \qquad
 dc_t=D_t\,dB_t-A_tc_t\,dt.
 \tag{2.4}
\]

Let

\[
 M_0=I,\qquad \dot M_t=A_tM_t.
 \tag{2.5}
\]

Since \(A_t\) is symmetric,

\[
 d(M_t^Tc_t)=M_t^TD_t\,dB_t.
 \tag{2.6}
\]

Thus \(M_t^Tc_t\) is a local martingale, and it is a true
square-integrable martingale in the compact/bounded regime of Section 1.
If

\[
 \alpha=\operatorname{Cov}_\mu(X,f)=c_0,
 \tag{2.7}
\]

then, in that regime,

\[
 \mathbb E[M_t^Tc_t]=\alpha,
 \qquad
 \mathbb E\langle M_tb,c_t\rangle=\langle b,\alpha\rangle.
 \tag{2.8}
\]

The orientation in (2.5) is essential: \(\dot M=A M\) gives
\(\dot M^T=M^TA\), which cancels the drift \(-Ac\) in (2.4).  No
commutativity of the matrices \(A_s\) is being assumed.

## 3. The exact consequence of (FA)

Assume now that \(\mu\) is isotropic, \(\mathbb Ef=0\), and
\(\mathbb Ef^2=1\).  Orthogonal projection onto linear functions gives

\[
 f=\alpha\cdot x+r,
 \qquad \mathbb E[Xr]=0,
 \qquad |\alpha|^2+\|r\|_2^2=1.
 \tag{3.1}
\]

Put \(\delta=\|r\|_2\) and, when \(\alpha\ne0\),
\(b=\alpha/|\alpha|\).  At time one, \(\mu_1\) is 1-strongly
log-concave.  Poincar\'e for \(\mu_1\), Cauchy--Schwarz for covariance,
and the posterior tower property give

\[
 |c_1|^2
 \le \mathbb E_1|\nabla f|^2,
 \qquad
 \mathbb E|c_1|^2
 \le \int|\nabla f|^2\,d\mu=:\mathcal E(f).
 \tag{3.2}
\]

Indeed, with \(u=c_1/|c_1|\),
\(|c_1|=\operatorname{Cov}_1(u\cdot X,f)\), while both
\(\operatorname{Var}_1(u\cdot X)\le1\) and
\(\operatorname{Var}_1(f)\le\mathbb E_1|\nabla f|^2\).

Suppose (FA) holds with constant \(C_{\mathrm{FA}}\):

\[
 \mathbb E\left[
   \frac{\langle M_1b,c_1\rangle^2}{|c_1|^2}
   \mathbf1_{\{c_1\ne0\}}
 \right]\le C_{\mathrm{FA}}.
 \tag{3.3}
\]

Using (2.8) and Cauchy--Schwarz gives the exact implication

\[
 \boxed{
 \mathcal E(f)\ge \frac{|\alpha|^2}{C_{\mathrm{FA}}}
 =\frac{1-\delta^2}{C_{\mathrm{FA}}}.}
 \tag{3.4}
\]

For a normalized first eigenfunction, \(\mathcal E(f)=\lambda\), so

\[
 \boxed{\lambda\ge (1-\delta^2)/C_{\mathrm{FA}}.}
 \tag{3.5}
\]

Without (FA), the same pairing and ordinary Cauchy--Schwarz yield

\[
 \mathbb E|M_1b|^2\ge\frac{1-\delta^2}{\lambda}.
 \tag{3.6}
\]

This corrects a possible overreading of the target: **(FA) by itself does
not give a universal lower bound on the whole spectral gap.**  It closes
only the near-linear branch, namely any branch on which
\(1-\delta^2=|\mathbb E[Xf]|^2\ge\eta>0\).  If \(\alpha=0\), (FA) has no
distinguished direction \(b\) and (3.5) gives no information.  A complete
KLS proof still needs a separate argument excluding genuinely nonlinear
small-gap modes.

The implication (3.4) does not require an attained eigenfunction; it
applies verbatim to every normalized locally Lipschitz test function for
which (2.8) is justified.  Hence it is compatible with a minimizing
sequence for the Poincar\'e quotient.

## 4. Endpoint law, Brownian bridges, and an exact static form of (FA)

The partition function satisfies

\[
 \partial_tZ=-\frac12\Delta_\theta Z.
 \tag{4.1}
\]

The transition density of (1.3), from \((s,h)\) to \((t,y)\), is the
Doob-transform kernel

\[
 p_{s,t}(h,y)=\varphi_{t-s}(y-h)\frac{Z(t,y)}{Z(s,h)},
 \tag{4.2}
\]

where \(\varphi_q\) is the centered Gaussian density with covariance
\(qI\).  In particular,

\[
 \rho_t(y)=\varphi_t(y)Z(t,y)
 \tag{4.3}
\]

is the density of \(\theta_t\), equivalently of
\(tX+\sqrt t\,G\).

In a finite-dimensional conditional path density, all factors \(Z\)
in (4.2) telescope.  Therefore

\[
 (\theta_s)_{0\le s\le t}\mid\{\theta_t=y\}
 \quad\hbox{is an ordinary Brownian bridge from \(0\) to \(y\).}
 \tag{4.4}
\]

This gives a precise endpoint formulation.  For unit \(b\), define

\[
 Q_{t,b}(y)=\mathbb E_{\mathrm{BB}(0,y)}\left[
   (M_tb)(M_tb)^T\right],
 \qquad
 M_t=\mathcal T\exp\!\left(\int_0^t
 A(s,\theta_s)\,ds\right).
 \tag{4.5}
\]

Then (FA) at time \(t=1\) is exactly

\[
 \boxed{
 \int \frac{\nabla F(1,y)^TQ_{1,b}(y)\nabla F(1,y)}
                {|\nabla F(1,y)|^2}
       \mathbf1_{\{|\nabla F(1,y)|>0\}}\rho_1(y)\,dy
 \le C.}
 \tag{4.6}
\]

Thus the fact that \(c_1\) is a posterior gradient does have exact
content: it selects a gradient direction in the Rayleigh quotient of a
Brownian-bridge-averaged product integral.  It does not, on its own,
control the top eigendirection of \(Q_{1,b}(y)\).

## 5. What endpoint transport determines, and what it does not

Let

\[
 v_{t,b}(y)=\mathbb E[M_tb\mid\theta_t=y],
 \qquad
 \Sigma_{t,b}(y)=Q_{t,b}(y)-v_{t,b}(y)v_{t,b}(y)^T\succeq0.
 \tag{5.1}
\]

Perturb the starting point of (1.3) from \(0\) to \(hb\).  Differentiating
(4.2) at \(h=0\), using that \(\mu\) is centered, gives for every smooth
compactly supported \(\phi\)

\[
 \int\langle\nabla\phi(y),v_{t,b}(y)\rangle\rho_t(y)\,dy
 =\frac1t\int \phi(y)\langle y,b\rangle\rho_t(y)\,dy.
 \tag{5.2}
\]

Equivalently,

\[
 -\operatorname{div}(\rho_tv_{t,b})
 =t^{-1}\langle y,b\rangle\rho_t
 \quad\hbox{in distributions}.
 \tag{5.3}
\]

Thus \(t v_{t,b}\) is a Stein-kernel column for the endpoint law.  This
identity recovers the mean pairing in (2.8), but it contains no information
about \(\Sigma_{t,b}\).  Indeed the (FA) integrand decomposes exactly as

\[
 \frac{c^TQ_{t,b}c}{|c|^2}
 =\frac{\langle c,v_{t,b}\rangle^2}{|c|^2}
  +\frac{c^T\Sigma_{t,b}c}{|c|^2}.
 \tag{5.4}
\]

The second term is nonnegative and is invisible to every endpoint
change-of-variables or first-variation formula.  This is a decisive
obstruction to deriving (FA) from coarea, the endpoint density, or the
Stein identity alone.  Such arguments determine a conditional mean; (FA)
requires a conditional second moment of a path-dependent, generally
noncommuting product integral.

## 6. The additional constraint supplied by a genuine eigenfunction

Assume now that \(d\mu=e^{-V}dx\) is regular enough for the weak
eigenvalue equation and that

\[
 -L_\mu f=\lambda f,
 \qquad \mathbb Ef=0,
 \qquad \mathbb Ef^2=1.
 \tag{6.1}
\]

For the observation \(Y_t=tX+\sqrt t\,G\), put

\[
 F_t(y)=\mathbb E[f(X)\mid Y_t=y],
 \qquad
 H_t(y)=\mathbb E[\nabla f(X)\mid Y_t=y].
 \tag{6.2}
\]

Testing the eigenvalue equation against
\(x\mapsto\mathbb E_G\phi(tx+\sqrt tG)\) gives

\[
 \lambda\int F_t\phi\,\rho_t
 =t\int\langle H_t,\nabla\phi\rangle\rho_t.
 \tag{6.3}
\]

Hence

\[
 -\operatorname{div}(\rho_tH_t)
 =\frac\lambda t\rho_tF_t,
 \qquad
 \int|H_t|^2\rho_t\le\int|\nabla f|^2d\mu=\lambda.
 \tag{6.4}
\]

When \(V\) and \(f\) are smooth, posterior integration by parts gives the
equivalent pointwise identity

\[
 \lambda F_t=(y-tm_t)\cdot H_t-t\operatorname{div}H_t,
 \qquad m_t(y)=\mathbb E[X\mid Y_t=y].
 \tag{6.5}
\]

Equations (4.6) and (6.4) are the exact formulation of the proposed
"posterior gradient of a fixed low mode" mechanism.  A proof of (FA) by
this route would have to establish the following bridge-alignment
principle, which is not implied by trace bounds:

> If \(F=\mathbb E[f(X)\mid X+G=\cdot]\) admits a flux \(H\) satisfying
> \(-\operatorname{div}(\rho H)=\lambda\rho F\) and
> \(\int|H|^2d\rho\le\lambda\), and if
> \(b=\mathbb E[Xf]/|\mathbb E[Xf]|\), then the Rayleigh quotient in
> (4.6) is universally bounded.

The missing step is specifically to control

\[
 \int
 \frac{\nabla F(y)^T\Sigma_{1,b}(y)\nabla F(y)}
      {|\nabla F(y)|^2}
 \mathbf1_{\{|\nabla F(y)|>0\}}\rho_1(y)\,dy.
 \tag{6.6}
\]

Neither (5.3) nor (6.4) contains \(\Sigma_{1,b}\).  Any claimed proof that
uses only those first-order divergence identities has dropped precisely
the path-fluctuation term (6.6).

## 7. It\^o audit of the normalized alignment functional

There is no hidden supermartingale behind the most direct normalization.
Let

\[
 v_t=M_tb,\qquad z_t=\langle v_t,c_t\rangle,
 \qquad r_t=|c_t|^2.
 \tag{7.1}
\]

From (2.4)--(2.6),

\[
 dz_t=\langle D_tv_t,dB_t\rangle,
 \qquad
 dr_t=2\langle D_tc_t,dB_t\rangle
 +\big(\|D_t\|_{\mathrm{HS}}^2-2c_t^TA_tc_t\big)dt.
 \tag{7.2}
\]

Stop before \(r_t\) reaches zero and before the quantities in (1.5)
diverge.  It\^o's formula gives the exact drift

\[
 \begin{split}
 d\left(\frac{z_t^2}{r_t}\right)_{\!\mathrm{drift}}
 ={}&\frac1{r_t}\left|D_t\left(v_t-
                    2\frac{z_t}{r_t}c_t\right)\right|^2dt\\
 &-\frac{z_t^2}{r_t^2}\|D_t\|_{\mathrm{HS}}^2dt
 +2\frac{z_t^2}{r_t^2}c_t^TA_tc_t\,dt.
 \end{split}
 \tag{7.3}
\]

The right side has no fixed sign.  In particular, when \(z_t=0\) but
\(D_tv_t\ne0\), its drift is strictly positive.  The first term contains
the already-amplified vector \(v_t\), so a trace estimate or an unweighted
bound on \(D_t\) does not close it.  A successful It\^o proof needs a new
compensator tied to the low-mode flux (6.4), or a genuinely new estimate
of the amplification-weighted Hessian term in (7.3).

## 8. Exact positive tests and the remaining model obstruction

The separate file `work/fa_counterexample_audit.md` proves (FA) without an
ambient trace estimate in the following cases:

1. isotropic Gaussian measure, where the quotient is exactly \(4\);
2. \(\operatorname{Unif}[-\sqrt3,\sqrt3]\otimes\gamma_{n-1}\) for its
   genuine first eigenfunction
   \(\sqrt2\sin(\pi x_1/(2\sqrt3))\), for which
   \(\delta^2=1-96/\pi^4<0.015\) and the quotient is at most \(e^6\)
   pathwise;
3. uniformly strongly log-concave measures;
4. arbitrary products of bounded-diameter blocks, with a constant
   independent of the number of blocks.

It also reduces the growing irreducible cone

\[
 p(s,z)\propto\mathbf1_{\{s>0\}}
 \exp\left(-s-\frac{|z|^2}{2s}\right)
 \tag{8.1}
\]

to an exact one-dimensional posterior and an arrowhead covariance matrix.
That model exhibits the remaining issue cleanly: the axial direction
couples to the instantaneous transverse tilt direction, which rotates
along the Brownian path, so the matrices \(A_t\) do not commute.  No
counterexample to (FA) was obtained, but fixed-endpoint posterior formulas
cannot settle it because they do not determine the bridge covariance
\(\Sigma_{1,b}\) in (6.6).

## 9. Verdict

The SDE and martingale ledger in the main note is correct in the regular
case, subject in the general case to an explicit uniform-integrability
passage.  The exact implication of (FA) is (3.5), not a full spectral-gap
bound unless a separate near-linearity dichotomy is supplied.

The posterior-gradient property produces the formalizable constraints
(4.6) and (6.4), but it does not remove the conditional path-variance term
(6.6).  This is the load-bearing unresolved estimate.  Trace bounds,
coarea, endpoint transport, and the normalized-ratio It\^o calculation do
not control it.  Any future proof should either:

* bound (6.6) using the low-mode flux \(H_t\) in a way that sees the same
  Brownian bridge as \(M_t\), or
* replace the synchronous product integral by a coupling whose conditional
  variance vanishes or is charged to the Bochner defect.

