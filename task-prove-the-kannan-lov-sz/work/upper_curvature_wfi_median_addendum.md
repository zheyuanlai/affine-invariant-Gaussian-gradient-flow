# Addendum: a self-contained proof of the one-dimensional (L^2) Stein bound

This addendum supplies the numerical estimate used in
`upper_curvature_wfi.md`:

\[
  \mathbb E\tau(S)^2\le 400
\]
for every centered, variance-one, one-dimensional log-concave law.  It is
intended to replace any appeal to an unproved ``standard'' median-density
constant.

Let (ho) be the density on its (possibly unbounded) interval (J),
let (m) be a median, and write (h=ho(m)).  Since the law is
continuous, both half-lines from (m) have mass (1/2).  Cantelli's
one-sided Chebyshev inequality gives (|m|le1): if (mge0), then
\[
 \frac12\le \mathbb P(S\ge m)le {1\over1+m^2};
\]
the case (mle0) follows by applying this to (-S).

Choose a mode (s_0).  If (s_0ge m), log-concavity makes (ho)
nondecreasing up to (s_0), so (ho(t)le h) for (tle m).  If
(s_0le m), the analogous statement holds on (tge m).  (When the
mode is an interval, choose either endpoint and use the corresponding
side.)  Thus one of the two half-lines has, after the change of variable
(u=|t-m|), a density (g) satisfying (0le gle h) and
(int_0^infty g(u),du=p:=1/2).

For completeness, the bathtub estimate used below is elementary.  Put
(L=p/h).  Since (u^2ge L^2) for (uge L),
\[
\begin{aligned}
 \int_0^infty u^2g(u),du
 &\ge \int_0^L u^2g(u),du+L^2\int_L^infty g(u),du\\
 &=L^2p-\int_0^L(L^2-u^2)g(u),du\\
 &\ge L^2p-h\int_0^L(L^2-u^2),du
  ={p^3\over3h^2}={1\over24h^2}.
\end{aligned}
\]
The selected side is part of the whole second moment about (m), hence
\[
  1+m^2=\mathbb E(S-m)^2\ge {1\over24h^2}.
\]
Together with (|m|le1), this gives (h\ge1/\sqrt{48}>1/8).

Let (T_+(s)=mathbb P(S\ge s)) and
(lambda_+(s)=ho(s)/T_+(s)).  The survival function of a
log-concave density is log-concave, so (lambda_+) is nondecreasing.
At the median (lambda_+(m)=2hge1/4).  Therefore, for (sge m),
(lambda_+(s)ge1/4), and
\[
 T_+(s+u)\le T_+(s)e^{-u/4},qquad u\ge0.
\]
Using the canonical Stein formula and
(ho(s)=lambda_+(s)T_+(s)),
\[
\begin{aligned}
 \tau(s)
 &= {sT_+(s)+\int_s^{\sup J}T_+(t)\,dt\over\rho(s)}\\
 &\le 4|s|+{4T_+(s)\over\rho(s)}
 \le 4|s|+16,qquad s\ge m.
\end{aligned}
\]
Here the first term is bounded by (4|s|) (if (s<0), it is
nonpositive and may simply be discarded), and the integral is at most
(4T_+(s)), while (T_+(s)/ho(s)le4).  The left-tail survival
(T_-(s)=mathbb P(Sle s)) has the analogous monotone reverse hazard;
the same calculation for (sle m) gives (	au(s)le4|s|+16).

Finally (mathbb E|S|le(mathbb ES^2)^{1/2}=1), so
\[
 \mathbb E\tau(S)^2
 \le \mathbb E(4|S|+16)^2
 \le 16+128+256=400.
\]
The argument is valid at hard endpoints by taking one-sided limits; all
quantities are nonnegative and the endpoint boundary terms in the Stein
formula vanish by the defining tail integrals.
