# Directional MMSE: scalar reduction, a weighted-Fisher gate, and an exact positive class

## 0. Verdict

Let \(X\) be centered, isotropic, and log-concave in \(\mathbb R^n\), let
\(G\sim N(0,I_n)\) be independent, and put

\[
 \varepsilon_b(X)=\mathbb E\operatorname {Var}(b\cdot X\mid X+G),
 \qquad |b|=1.
 \tag{0.1}
\]

No unconditional proof of

\[
 \inf_{n,X,b}\varepsilon_b(X)>0
 \tag{DMMSE}
\]

and no log-concave counterexample is obtained here.  The logical obstruction
is exact: KLS implies (DMMSE), while a counterexample to (DMMSE) would refute
KLS.  The new concrete outputs of this audit are the following.

1. There is a fully rigorous scalar reduction, including conditional
   log-concavity, which gives

   \[
   \varepsilon_b(X)
   \ge {3\over44800}
       {d\over \log _2^2(128/d)},
   \qquad
   d=\mathbb E\operatorname {Var}(b\cdot X\mid
                    P_{b^\perp}X+G_{b^\perp}).
   \tag{0.2}
   \]

2. If \(C=C_P(P_{b^\perp}X+G_{b^\perp})\), then

   \[
   d\ge{1\over1+C},
   \quad
   \varepsilon_b(X)\ge
   {3\over44800(1+C)\log _2^2(128(1+C))}.
   \tag{0.3}
   \]

   This proves the proposed dimension-descent estimate, with explicit
   constants, but it is not dimension free unless a new estimate for the
   \((n-1)\)-dimensional marginal is supplied.

3. A sharper sufficient condition is a one-dimensional weighted Fisher
   bound for the noisy orthogonal slice channel.  If \(S=b\cdot X\),
   \(Z=P_{b^\perp}X+G_{b^\perp}\), \(\tau\) is the canonical Stein kernel of
   \(S\), and

   \[
   \mathcal I_\perp
   :=\mathbb E\left[\tau(S)^2
       \left|\partial_s\log q_s(Z)\big|_{s=S}\right|^2\right]\le K,
   \tag{0.4}
   \]

   where \(q_s\) is the density of \(Z\mid S=s\), then

   \[
   \varepsilon_b(X)\ge {1\over4(400+K)}.
   \tag{0.5}
   \]

4. The weighted Fisher bound holds with \(K=1\), in every dimension, for
   the complete zero-mean conditional-Gaussian class

   \[
   p(t,y)\propto
   \exp\left[-W(t)-\frac12y^TQ(t)y\right],
   \qquad R(t)=Q(t)^{-1},
   \tag{0.6}
   \]

   whenever the joint density is log-concave, \(\operatorname {Var}T=1\),
   and \(\mathbb E R(T)=I\).  Noncommuting and rotating matrices \(R(t)\)
   are allowed.  Consequently this class satisfies

   \[
   \varepsilon_{e_T}(X)\ge {1\over1604}
   \tag{0.7}
   \]

5. A smooth isotropic non-log-concave family has
   \(\varepsilon_{e_1}\to0\).  It keeps a linear trace of posterior error
   and a nonvanishing entropy power.  Thus scalar marginal MMSE, isotropic
   covariance, and entropy/trace estimates cannot by themselves imply the
   directional statement.  The missing ingredient must use joint
   log-concavity to prohibit nonlinear coding in the orthogonal coordinates.

The remaining general target is therefore the weighted Fisher inequality
(0.4), or a substitute that controls precisely the same nonlinear slice
information.  It is verified below for (0.6), but no argument extending it
to arbitrary log-concave slices is known.  In differential form that
extension is a reinforced Prekopa/Brascamp--Lieb stability statement, not a
scalar-channel estimate.

---

## 1. Exact channel identities and the KLS implication

Write \(Y=X+G\), let \(q\) be its density, and set

\[
 m(y)=\mathbb E[X\mid Y=y],\qquad
 A(y)=\operatorname {Cov}(X\mid Y=y).
 \tag{1.1}
\]

The posterior has potential

\[
 V(x)+\frac12|x-y|^2,
\]

so it is 1-strongly log-concave on its convex support.  Brascamp--Lieb and
exponential-family differentiation give

\[
 0\preceq A(y)\preceq I,
 \qquad Dm(y)=A(y),
 \qquad \nabla\log q(y)=m(y)-y.
 \tag{1.2}
\]

For \(h_b(y)=b\cdot m(y)\), total variance gives

\[
 \operatorname {Var}_q h_b=1-\varepsilon_b,
 \qquad
 \int|\nabla h_b|^2dq
 =\mathbb E|A(Y)b|^2\le\varepsilon_b.
 \tag{1.3}
\]

Poincare for \(q\) therefore implies

\[
 \varepsilon_b\ge {1\over1+C_P(q)}.
 \tag{1.4}
\]

The convolution inequality \(C_P(\mu*\gamma)\le C_P(\mu)+1\) follows by
conditioning a test function first on \(G\), then on \(X\), applying the
two Poincare inequalities, and using Jensen for the conditional gradients.
Hence

\[
 \boxed{\quad
 \varepsilon_b\ge {1\over C_P(\mu)+2}.
 \quad}
 \tag{1.5}
\]

Thus KLS proves (DMMSE).  Conversely, if a sequence has
\(\varepsilon_{b_n}\to0\), (1.3) gives

\[
 C_P(\mu_n*\gamma)\ge{1-\varepsilon_{b_n}\over\varepsilon_{b_n}}
 \longrightarrow\infty,
 \tag{1.6}
\]

and the convolution inequality forces \(C_P(\mu_n)\to\infty\).  Therefore
a log-concave counterexample to (DMMSE) is automatically a counterexample
to KLS.

In Fisher notation, integration by parts gives

\[
 \mathbb EA(Y)=I-J(q),
 \qquad
 \varepsilon_b=1-b^TJ(q)b.
 \tag{1.7}
\]

An operator upper bound \(J(q)\preceq(1-c)I\) is exactly (DMMSE), not an
independent generic Fisher inequality.

---

## 2. The conditional scalar reduction

Rotate so that \(b=e_1\), and write

\[
 S=X_1,qquad X_\perp=P_{e_1^\perp}X,qquad
 Z=X_\perp+G_\perp,qquad U=S+G_1.
 \tag{2.1}
\]

The full observation is \((U,Z)\).  Put

\[
 h(z)=\mathbb E[S\mid Z=z],qquad
 v(z)=\operatorname {Var}(S\mid Z=z),qquad
 d=\mathbb Ev(Z).
 \tag{2.2}
\]

### Lemma 2.1: conditional log-concavity

For \(q_Z\)-almost every \(z\), the conditional law
\(\nu_z=\mathcal L(S\mid Z=z)\) is log-concave on \(\mathbb R\).

**Proof.**  On the affine hull of \(\mu\), write its density as
\(e^{-V}\) with \(V\) convex.  The joint density of \((S,Z)\) is

\[
 r(s,z)=\int_{e_1^\perp}
 e^{-V(se_1+x)}\,\gamma_{n-1}(z-x)\,dx.
 \tag{2.3}
\]

The integrand is jointly log-concave in \((s,z,x)\).  Prekopa's theorem
therefore makes \(r\) log-concave in \((s,z)\).  Dividing the one-dimensional
slice \(r(\cdot,z)\) by its integral preserves log-concavity.  The argument
is intrinsic if the original support is an affine subspace.  Since an
isotropic law has full-dimensional affine hull, no ambient degeneracy
remains in the present setting. \(\square\)

### Lemma 2.2: scalar unit-noise MMSE

If \(T\) is a one-dimensional log-concave random variable of variance
\(v\), and \(N\sim N(0,1)\) is independent, then

\[
 \mathbb E\operatorname {Var}(T\mid T+N)
 \ge {v\over12v+2}
 \ge {1\over14}\min(v,1).
 \tag{2.4}
\]

**Proof.**  Let \(Q=T+N\), \(a(Q)=\mathbb E[T\mid Q]\), and denote the
MMSE by \(e\).  The one-dimensional log-concave Poincare bound is
\(C_P(T)\le12v\), and convolution gives \(C_P(Q)\le12v+1\).  The posterior
of \(T\) given \(Q=q\) is 1-strongly log-concave, so

\[
 a'(q)=\operatorname {Var}(T\mid Q=q)\in[0,1].
 \tag{2.5}
\]

Thus

\[
 \operatorname {Var}a(Q)=v-e,
 \qquad
 \mathbb E a'(Q)^2\le\mathbb E a'(Q)=e.
 \tag{2.6}
\]

Poincare for \(Q\) yields \(v-e\le(12v+1)e\), which is (2.4).
The argument extends to interval supports by approximation. \(\square\)

Conditioning first on \(Z\), Lemmas 2.1--2.2 give the exact reduction

\[
 \boxed{\quad
 \varepsilon_{e_1}(X)
 \ge {1\over14}\mathbb E\min(v(Z),1).
 \quad}
 \tag{2.7}
\]

This is the precise point at which a scalar-only argument stops: it still
has to show that noisy orthogonal coordinates do not make \(v(Z)\) tiny.

### Lemma 2.3: a uniform tail conversion

For \(S\) centered, variance one, and log-concave, the standard
one-dimensional Borell moment estimate may be normalized as

\[
 \|S\|_r\le2r,\qquad r\ge2.
 \tag{2.8}
\]

Consequently, for every \(p\ge1\),

\[
 \|v(Z)\|_p
 \le\|\mathbb E[S^2\mid Z]\|_p
 \le\|S^2\|_p
 =\|S\|_{2p}^2
 \le16p^2.
 \tag{2.9}
\]

Let \(0<d=\mathbb Ev\le1\), put

\[
 L=\log _2(128/d),\qquad p=\lceil4L\rceil,
 \qquad R=32p^2.
 \tag{2.10}
\]

By (2.9),

\[
 \mathbb E[v;v>R]
 \le {\mathbb Ev^p\over R^{p-1}}
 \le R,2^{-p}.
 \tag{2.11}
\]

Since \(L\ge7\), \(p\le4L+1\le5L\), and
\((4L+1)^2\le2^{3L}\).  Hence

\[
 R2^{-p}\le32\,2^{3L}2^{-4L}=32\,2^{-L}=d/4.
 \tag{2.12}
\]

As \(R\ge1\),

\[
 \mathbb E\min(v,1)
 \ge {1\over R}\mathbb E[v;v\le R]
 \ge {3d\over4R}
 \ge {3d\over3200L^2}.
 \tag{2.13}
\]

Combining (2.7) and (2.13) proves (0.2).

---

## 3. The exact dimension-descent estimate

The marginal density of \(X_\perp\) is log-concave.  Given \(Z=z\), the
marginal posterior of \(X_\perp\) has density proportional to

\[
 p_\perp(x)\exp(-|x-z|^2/2),
 \tag{3.1}
\]

and is therefore 1-strongly log-concave.  Its covariance matrix \(C(z)\)
satisfies \(C(z)\preceq I\).

Exponential-family differentiation gives

\[
 \nabla h(z)=\operatorname {Cov}(S,X_\perp\mid Z=z)=:c(z).
 \tag{3.2}
\]

For \(c(z)\ne0\), conditional Cauchy--Schwarz and \(C(z)\preceq I\) give

\[
 |c|^4
 =\operatorname {Cov}(S,c\cdot X_\perp\mid Z)^2
 \le v\operatorname {Var}(c\cdot X_\perp\mid Z)
 \le v|c|^2.
 \tag{3.3}
\]

Thus

\[
 |\nabla h(z)|^2\le v(z).
 \tag{3.4}
\]

The law of total variance and Poincare for \(q_Z\) now yield

\[
 1=d+\operatorname {Var}_{q_Z}h
 \le d+C_P(q_Z)\mathbb E|\nabla h|^2
 \le(1+C_P(q_Z))d.
 \tag{3.5}
\]

This proves \(d\ge(1+C_P(q_Z))^{-1}\).  Substitution in (0.2), noting that
\(d\mapsto d/\log_2^2(128/d)\) is increasing, proves (0.3).

If \(K_{n-1}\) denotes the supremum of the Poincare constants of isotropic
log-concave laws in dimension \(n-1\), then

\[
 C_P(q_Z)\le C_P(PX)+1\le K_{n-1}+1.
 \tag{3.6}
\]

Therefore

\[
 \varepsilon_b(X)
 \ge {3\over44800(K_{n-1}+2)
             \log_2^2(128(K_{n-1}+2))}.
 \tag{3.7}
\]

Equation (3.7) is a proved recursion, not a dimension-free conclusion.  A
dimension-free insertion for \(K_{n-1}\) is KLS itself.

---

## 4. A weighted Fisher-information gate

Let \(S\) be any centered, variance-one, one-dimensional log-concave
variable with density \(\rho=e^{-\Phi}\) on its support interval.  Its
canonical Stein kernel is

\[
 \tau(s)={1\over\rho(s)}\int_s^\infty u\rho(u)\,du.
 \tag{4.1}
\]

It satisfies

\[
 \mathbb E[S\varphi(S)]=\mathbb E[\tau(S)\varphi'(S)],
 \qquad \mathbb E\tau(S)=1.
 \tag{4.2}
\]

For a regular channel \(W\mid S=s\) with density \(r_s(w)\), write

\[
 \ell_s(w)=\partial_s\log r_s(w),
 \qquad
 I_W(s)=\int\ell_s(w)^2r_s(w)\,dw.
 \tag{4.3}
\]

### Lemma 4.1: Stein--Fisher MMSE inequality

Let \(a(W)=\mathbb E[S\mid W]\) and
\(D=\mathbb E(S-a(W))^2\).  Then

\[
 \boxed{\quad
 (1-D)^2\le D\,\mathbb E[\tau(S)^2I_W(S)].
 \quad}
 \tag{4.4}
\]

**Proof.**  Orthogonality of conditional expectation gives

\[
 D=\mathbb E[S(S-a(W))].
 \tag{4.5}
\]

For fixed \(w\), apply (4.2) to
\(s\mapsto r_s(w)(s-a(w))\), then integrate in \(w\).  Since
\(\partial_sr_s=r_s\ell_s\),

\[
 D=\mathbb E\tau(S)
   +\mathbb E[\tau(S)(S-a(W))\ell_S(W)]
 =1+\mathbb E[\tau(S)(S-a(W))\ell_S(W)].
 \tag{4.6}
\]

Cauchy--Schwarz proves (4.4). \(\square\)

### Lemma 4.2: explicit \(L^2\) bound for the Stein kernel

For a centered variance-one log-concave density,

\[
 \mathbb E\tau(S)^2\le400.
 \tag{4.7}
\]

**Proof.**  A standard sharp one-dimensional log-concavity fact is that if
\(m\) is a median, then \(\rho(m)\ge1/8\) when the variance is one.  Also
\(|m|\le\sqrt2\) by Chebyshev.  The right and left hazard rates of a
log-concave density are monotone.  At the median each is at least \(1/4\).
For \(s\ge m\), integration of the survival function therefore gives

\[
 \tau(s)
 ={s\,\mathbb P(S\ge s)+\int_s^\infty\mathbb P(S\ge u)du\over\rho(s)}
 \le4|s|+16.
 \tag{4.8}
\]

The analogous left-tail formula gives the same bound for \(s\le m\).
Consequently

\[
 \mathbb E\tau^2
 \le\mathbb E(4|S|+16)^2
 \le16+128\mathbb E|S|+256\le400.
 \tag{4.9}
\]

The median-density constant \(1/8\) is deliberately nonsharp; it follows,
for example, from the elementary one-dimensional bounds relating median
density and standard deviation. \(\square\)

Return now to \(S=X_1\) and \(Z=X_\perp+G_\perp\).  Let \(q_s\) be the
conditional density of \(Z\mid S=s\), and define

\[
 K=\int\rho(s)\tau(s)^2
      \int|\partial_s\log q_s(z)|^2q_s(z)\,dz\,ds.
 \tag{4.10}
\]

For the full observation \(W=(S+G_1,Z)\), the conditional score is

\[
 \ell_s(u,z)=(u-s)+\partial_s\log q_s(z).
 \tag{4.11}
\]

The two summands are conditionally independent and centered, so

\[
 \mathbb E[\tau(S)^2I_W(S)]
 =\mathbb E\tau(S)^2+K\le400+K.
 \tag{4.12}
\]

If \(D\le1/2\), (4.4) gives \(1/4\le(400+K)D\); if \(D>1/2\) the same
lower bound is automatic.  Therefore

\[
 \boxed{\quad
 \varepsilon_{e_1}(X)=D\ge {1\over4(400+K)}.
 \quad}
 \tag{4.13}
\]

This proves (0.5).

The isotropy cancellation relevant to (4.10) is exact.  If
\(m_\perp(s)=\mathbb E[X_\perp\mid S=s]\), then

\[
 \mathbb E[\tau(S)\partial_s\log q_S(Z)]=0,
 \tag{4.14}
\]

and

\[
 \mathbb E[\tau(S)Z\,\partial_s\log q_S(Z)]
 =\mathbb E[\tau(S)m_\perp'(S)]
 =\mathbb E[SX_\perp]=0.
 \tag{4.15}
\]

Thus the weighted slice score is orthogonal to constants and to every
linear transverse output coordinate.  Bounding its remaining nonlinear
\(L^2\) mass is exactly the live issue.

---

## 5. Complete proof for noncommuting conditional-Gaussian slices

Consider the density

\[
 p(t,y)=Z^{-1}\exp\left[-W(t)-\frac12y^TQ(t)y\right],
 \qquad R(t)=Q(t)^{-1}\succ0,
 \tag{5.1}
\]

on an interval times \(\mathbb R^m\).  Assume it is centered and isotropic:

\[
 \mathbb ET=0,qquad \operatorname {Var}T=1,qquad
 \mathbb E R(T)=I_m.
 \tag{5.2}
\]

There are no cross covariances because \(\mathbb E[Y\mid T]=0\).

The block Hessian Schur complement shows that (5.1) is jointly
log-concave if and only if

\[
 W''(t)\ge0,
 \qquad R''(t)\preceq0.
 \tag{5.3}
\]

Indeed,

\[
 Q''-2Q'Q^{-1}Q'=-QR''Q.
 \tag{5.4}
\]

The marginal density of \(T\) is proportional to \(e^{-\Phi(t)}\), with

\[
 \Phi(t)=W(t)-\frac12\log\det R(t)+\text{constant}.
 \tag{5.5}
\]

Let \(Z=Y+G_m\).  Given \(T=t\),

\[
 Z\sim N(0,R(t)+I).
 \tag{5.6}
\]

The Fisher information of this conditional experiment with respect to
the scalar parameter \(t\) is

\[
 I_Z(t)=\frac12
 \left\|(R+I)^{-1/2}R'(R+I)^{-1/2}\right\|_{\mathrm{HS}}^2.
 \tag{5.7}
\]

Since \(R^{1/2}(R+I)^{-1/2}\) is a contraction, even when \(R\) and
\(R'\) do not commute,

\[
 I_Z(t)\le\frac12
 \left\|R^{-1/2}R'R^{-1/2}\right\|_{\mathrm{HS}}^2.
 \tag{5.8}
\]

On the other hand, (5.3) gives

\[
\begin{aligned}
 \Phi''(t)
 &=W''(t)+\frac12\operatorname {Tr}
   (R^{-1}R'R^{-1}R'-R^{-1}R'')\\
 &\ge\frac12
 \left\|R^{-1/2}R'R^{-1/2}\right\|_{\mathrm{HS}}^2.
\end{aligned}
 \tag{5.9}
\]

Hence

\[
 \boxed{\quad I_Z(t)\le\Phi''(t).\quad}
 \tag{5.10}
\]

There is an exact Stein-curvature identity.  From
\((\tau\rho)'=-t\rho\),

\[
 \tau'(t)-\tau(t)\Phi'(t)=-t.
 \tag{5.11}
\]

Using \(\mathbb Eg'=\mathbb E[g\Phi']\) with
\(g=\tau^2\Phi'\), one obtains

\[
\begin{aligned}
 \mathbb E[\tau^2\Phi'']
 &=\mathbb E[\tau^2(\Phi')^2-2\tau\tau'\Phi']\\
 &=\mathbb E[T^2-(\tau')^2]
 =1-\mathbb E(\tau')^2\le1.
\end{aligned}
 \tag{5.12}
\]

All boundary terms vanish first for smooth confined densities; monotone
approximation gives the interval and extended-valued convex cases.  From
(5.10)--(5.12),

\[
 K=\mathbb E[\tau(T)^2I_Z(T)]\le1.
 \tag{5.13}
\]

Equation (4.13) now proves

\[
 \boxed{\quad
 \mathbb E\operatorname {Var}(T\mid T+G_1,Y+G_m)
 \ge {1\over1604}.
 \quad}
 \tag{5.14}
\]

No simultaneous diagonalization of the matrices \(R(t)\) is used.  Thus
rotating, noncommuting covariance eigenspaces are fully included.

---

## 6. A smooth isotropic coding family without log-concavity

This construction proves that the joint log-concavity hypothesis, rather
than scalar log-concavity or covariance alone, must prevent nonlinear
orthogonal coding.

Fix \(0<\delta<1\).  Let \(B\) be uniform on \(\{-1,1\}\).  Conditional on
\(B=b\), let \(X_1,\ldots,X_m\) be independent centered Gaussians of
variance \(1+b\delta\).  Let \(H\sim N(0,1)\) be independent and put

\[
 X_0={B+\eta H\over\sqrt{1+\eta^2}}.
 \tag{6.1}
\]

Then \(X^{(m,\eta)}=(X_0,X_1,\ldots,X_m)\) has a smooth positive density
and covariance \(I_{m+1}\): all means vanish, every coordinate has variance
one, and every cross covariance is zero.  The density is a two-component
Gaussian mixture and is not log-concave for the parameter regime used
below.

Observe \(Y_i=X_i+G_i\).  Given \(B=b\), the transverse outputs are iid
\(N(0,2+b\delta)\).  The classifier

\[
 \widehat B=operatorname {sign}\left(rac1m\sum_{i=1}^mY_i^2-2\right)
 \tag{6.2}
\]

has, by Chebyshev,

\[
 \mathbb P(\widehat B\ne B\mid B=b)
 \le {2(2+b\delta)^2\over m\delta^2}
 \le {2(2+\delta)^2\over m\delta^2}.
 \tag{6.3}
\]

Consequently

\[
 \mathbb E\operatorname {Var}(B\mid Y_1,\ldots,Y_m)
 \le\mathbb E(B-\widehat B)^2
 \le {8(2+\delta)^2\over m\delta^2}.
 \tag{6.4}
\]

The estimator \(\widehat X_0=\widehat B/\sqrt{1+\eta^2}\) therefore gives

\[
\begin{aligned}
 \mathbb E\operatorname {Var}(X_0\mid Y_0,ldots,Y_m)
 &\le\mathbb E(X_0-\widehat X_0)^2\\
 &\le {2\over1+\eta^2}
   \left(4\mathbb P(\widehat B\ne B)+\eta^2\right)\\
 &\le {64(2+\delta)^2\over m\delta^2}+2\eta^2.
\end{aligned}
 \tag{6.5}
\]

Taking \(\eta=m^{-1/2}\) makes the directional MMSE tend to zero.  The
direct observation \(Y_0\) can only reduce it further.

At the same time, conditioning additionally on \(B\) can only decrease
posterior variance, and then the transverse posterior factorizes.  Hence

\[
 \mathbb E\operatorname {Var}(X_i\mid Y_0,\ldots,Y_m)
 \ge\mathbb E\operatorname {Var}(X_i\mid Y_0,\ldots,Y_m,B)
 =\mathbb E{1+B\delta\over2+B\delta}
 \ge {1-\delta\over2-\delta}.
 \tag{6.6}
\]

Thus \(\operatorname {Tr}\mathbb E A(Y)\asymp m\).  Moreover,
the loss in the input differential entropy caused by
\(\eta=m^{-1/2}\) is only \(O(\log m)\); divided by \(m+1\), it vanishes.
Thus the entropy power stays bounded below by a numerical constant while
the least directional posterior error vanishes.  This is a realized,
smooth version of the algebraic trace-to-operator obstruction.

---

## 7. Exact gap for arbitrary log-concave slices

For the general isotropic log-concave law, define

\[
 q_s(z)=\text{density of }X_\perp+G_\perp\mid S=s,
 \qquad
 I_\perp(s)=\int|\partial_s\log q_s|^2q_s.
 \tag{7.1}
\]

The narrow sufficient statement exposed by Section 4 is

\[
 \boxed{\qquad
 \int\tau(s)^2I_\perp(s)\rho(s)\,ds\le C.
 \qquad}
 \tag{WFI}
\]

Its hypotheses include the exact isotropy cancellations (4.14)--(4.15).
If (WFI) holds, (4.13) proves (DMMSE) with
\(c=1/[4(400+C)]\).

A tempting pointwise inequality

\[
 I_\perp(s)\le\Phi''(s)
 \tag{7.2}
\]

is false for arbitrary joint log-concave families.  For the affine Gaussian
location model \(X_\perp=aS+\sigma G\), the left side after unit Gaussian
noise is \(|a|^2/(1+\sigma^2)\), while a standard Gaussian prior has
\(\Phi''=1\); \(|a|\) is arbitrary.  This example has
\(\operatorname {Cov}(S,X_\perp)=a\ne0\), so isotropy removes exactly this
affine obstruction.  What remains is a nonlinear stability problem.

For a smooth joint potential \(V(s,x)\), Prekopa gives

\[
 \Phi''(s)
 =\mathbb E_sV_{ss}-\operatorname {Var}_s(V_s)\ge0.
 \tag{7.3}
\]

Convexity and Brascamp--Lieb prove the nonnegativity through

\[
 \operatorname {Var}_s(V_s)
 \le\mathbb E_s\langle V_{xx}^{-1}V_{xs},V_{xs}\rangle
 \le\mathbb E_sV_{ss}.
 \tag{7.4}
\]

To establish (WFI), one needs a stability version of (7.4) after the
affine score component has been removed by (4.15), and in the metric seen
after unit Gaussian convolution.  Equality in (7.4) contains affine
Gaussian translations; isotropy eliminates those, but no proved general
remainder controls all nonlinear slice motion.  This is the same
affine-Hessian rigidity obstruction that appears in the low-mode
mean-gradient route.

The audits of the standard alternatives are therefore:

* **Scalar MMSE:** proves (2.7), but cannot control \(v(Z)\).
* **Poincare of the transverse output:** proves (0.3), but a uniform input
  is KLS in dimension \(n-1\).
* **Entropy power:** controls \(\operatorname {Tr}\mathbb EA\), not its
  least eigenvalue; Section 6 realizes the gap.
* **Fisher information:** \(J(\mu*\gamma)\preceq I\) is automatic, but a
  uniform strict operator gap is exactly (DMMSE).
* **Pointwise posterior covariance:** false even in one dimension near a
  hard endpoint; only the averaged target is viable.
* **Unsmoothened fiber-width bounds:** false for isotropic cubes in a
  diagonal direction, where the projection-volume/directional-variation
  ratio is of order \(\sqrt n\).  Unit transverse Gaussian smoothing is
  essential.

No step above assumes KLS or a dimension-free Poincare inequality.  The
only unproved general statement is explicitly isolated as (WFI).
