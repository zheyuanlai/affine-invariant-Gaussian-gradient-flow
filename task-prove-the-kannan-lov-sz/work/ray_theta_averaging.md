# Averaging the Cheeger--Jensen direction on calibrated rays

## 1. A quantified boundary--quotient lemma

Let \(S\) be a half-mass finite-perimeter cut and suppose that, on its
transport region, there is a calibrated-ray disintegration

\[
 d\mu(x)=d\nu_y(t)\,d\eta(y),\qquad
 x=x_y(t),\qquad S=\{t>0\}.                           \tag{1.1}
\]

Assume every \(\nu_y\) is one-dimensional log-concave and balanced at zero:
\(\nu_y(t>0)=\nu_y(t<0)=1/2\). Write

\[
 b_y={d\nu_y\over dt}(0),\qquad
 p=P_\mu(S)=\int b_y\,d\eta(y),                       \tag{1.2}
\]

and let \(N_y=\dot x_y(t)\) be the oriented ray direction. These hypotheses
hold for the balanced calibrated quotient of a threshold Kantorovich
potential; the statement below is conditional on the Cheeger cut and that
threshold quotient referring to the same set.

For the Gaussian channel at variance \(s\), retain the notation

\[
 g_s(Y)=\mathbb E[\mathbf1_S(X)\mid Y],\quad
 W(Y)=\nabla g_s(Y),\quad
 h_s(X)=T_s\mathbf1_S(X)=\mathbb E[g_s(Y)\mid X].
\]

Put

\[
 m_s(X)=\nabla h_s(X)=\mathbb E[W(Y)\mid X],\qquad
 \theta_s(X)={m_s(X)\over|m_s(X)|}                    \tag{1.3}
\]

where the last vector is arbitrary when \(m_s=0\), and define the raywise
approximation error

\[
 E_y(s)=\int|h_s(x_y(t))-\mathbf1_{\{t>0\}}|\,d\nu_y(t).
                                                               \tag{1.4}
\]

The global error is exact:

\[
             \int E_y(s)d\eta(y)=\int|h_s-\mathbf1_S|d\mu=2U(s).
                                                               \tag{1.5}
\]

The mixed boundary--quotient error is

\[
                    \mathfrak B_s=\int b_yE_y(s)d\eta(y).       \tag{1.6}
\]

**Lemma 1.1.** Under (1.1)--(1.6),

\[
 \boxed{
 \int |m_s(X)|\left(1-
       |\langle\theta_s(X),N_{Q(X)}\rangle|\right)d\mu(X)
 \le2\mathfrak B_s.}                                  \tag{1.7}
\]

Equivalently, for the unoriented projectors,

\[
 \boxed{
 \int |m_s(X)|
  \|\theta_s\theta_s^T-N_QN_Q^T\|_{HS}^2d\mu(X)
 \le8\mathfrak B_s.}                                  \tag{1.8}
\]

**Proof.** The one-dimensional Cheeger constant of a log-concave law
balanced at zero is \(2b_y\). Its \(L^1\) Cheeger inequality and the triangle
inequality give

\[
 \begin{aligned}
 \int\left|{d\over dt}h_s(x_y(t))\right|d\nu_y(t)
 &\ge2b_y\inf_c\int|h_s(x_y(t))-c|d\nu_y(t)\\
 &\ge b_y\{1-2E_y(s)\}.                               \tag{1.9}
 \end{aligned}
\]

On a calibrated ray,

\[
 {d\over dt}h_s(x_y(t))=\langle m_s(x_y(t)),N_y\rangle.
\]

Integrating (1.9) over \(y\), using (1.2), and using the general heat-profile
upper bound

\[
                    \int|m_s|d\mu\le H(s)\le p
\]

gives

\[
 \int\{|m_s|-|\langle m_s,N_Q\rangle|\}d\mu
 \le p-\{p-2\mathfrak B_s\}=2\mathfrak B_s,
\]

which is (1.7). Finally,

\[
 \|aa^T-bb^T\|_{HS}^2
 =2(1-\langle a,b\rangle^2)
 \le4(1-|\langle a,b\rangle|)
\]

for unit \(a,b\), proving (1.8). \(\square\)

Define the amplitude-weighted ray average

\[
 A_y=\int|m_s(x_y(t))|d\nu_y(t),\qquad
 \overline P_y={1\over A_y}\int |m_s|\,
                    \theta_s\theta_s^T\,d\nu_y               \tag{1.10}
\]

when \(A_y>0\). Convexity of the squared Hilbert--Schmidt norm and (1.8)
give

\[
 \boxed{
 \int A_y\|\overline P_y-N_yN_y^T\|_{HS}^2d\eta(y)
 \le8\mathfrak B_s.}                                  \tag{1.11}
\]

Thus averaging \(\theta\) along the rays does recover the calibrated normal
line, but only with the mixed weight \(\mathfrak B_s\).

## 2. Combination with the audited Jensen direction

Let \(u(Y)=W(Y)/|W(Y)|\), and let

\[
 \delta_J(s)=\mathbb E|W|-\mathbb E|\mathbb E[W\mid X]|.
\]

Section 12 of heatflow_bernstein.md proves

\[
 \mathbb E[|W||u-\theta_s(X)|^2]=2\delta_J(s).        \tag{2.1}
\]

Combining (2.1) with (1.8) also aligns the posterior active line directly
with the calibrated ray line:

\[
 \boxed{
 \mathbb E\!\left[
 |W(Y)|\|u(Y)u(Y)^T-N_QN_Q^T\|_{HS}^2\right]
 \le12\delta_J(s)+16\mathfrak B_s.}                  \tag{2.2}
\]

To verify the constants, first use
\(\|uu^T-\theta\theta^T\|_{HS}^2\le2|u-\theta|^2\).
Next put \(a(X)=\mathbb E[|W|\mid X]\). The passage from the
\(|m_s|\)-weight in (1.8) to the \(a\)-weight costs at most
\(2\int(a-|m_s|)=2\delta_J\), because a projector distance squared is at
most two. Finally apply the squared triangle inequality for the two
projector differences. This gives \(8\delta_J+4\delta_J+16\mathfrak B_s\).

For a Cheeger \(\varepsilon\)-minimizer, the audited estimate is

\[
             \delta_J(s)\le\varepsilon+4pU(s).        \tag{2.3}
\]

If all ray laws have standard deviation at least \(L\), the standard
one-dimensional log-concave density bound
\(\|d\nu_y/dt\|_\infty\le1/\sigma_y\) gives

\[
 \mathfrak B_s\le {1\over L}\int E_y\,d\eta
 ={2U(s)\over L}.                                    \tag{2.4}
\]

Equations (2.2)--(2.4) are a genuine dimension-free line-calibration lemma.
They do not localize to a selected long-ray band: on the complement,
\(b_yE_y\) can dominate (1.6), and global Cheeger minimality supplies no
packetwise upper bound on \(\int|m_s|\). This boundary--quotient correlation
is the remaining \(B\)--\(Q\) synergy term.

Nor does (1.11) create a quotient height for the finite Bregman competitor.
It controls an amplitude-weighted average of projectors, but gives no
Lipschitz or switching-cost control for a scalar function of \(Q\).

## 3. Sharp log-concave countermodel for the amplitude term

The amplitude term cannot be controlled by ray averaging, even when all
directions agree exactly. Let

\[
 \mu_K=N(0,K)\quad\hbox{on }\mathbb R,\qquad
 S=(0,\infty),\qquad 0<\alpha\le {1\over4},\qquad s=\alpha K.
                                                               \tag{3.1}
\]

This is an exact balanced Cheeger halfspace. Its calibrated-ray quotient has
one infinite ray, with direction \(N=+1\) and conditional scale \(\sqrt K\).
For every \(y\),

\[
 g_s(y)=\Phi\!\left(
 {y\over\sqrt{\alpha(1+\alpha)K}}\right),\qquad W(y)=g_s'(y)>0.
                                                               \tag{3.2}
\]

Consequently

\[
                 u(Y)=\theta_s(X)=N=+1                 \tag{3.3}
\]

identically. The Jensen angular deficit, the ray-averaging deficit, and the
phase term in (12.24) are all exactly zero.

Nevertheless the amplitude term at the natural regularization
\(\lambda=M_s^2\), \(M_s=I(1/2)/\sqrt s\), has the lower bound

\[
 \boxed{
 D_{\rm amp}:={1\over2}\mathbb E[
   (\kappa_\lambda(|W(Y)|)-\kappa_\lambda(|W(Y')|))^2]
 \ge 5\cdot10^{-6}\sqrt\alpha,}                     \tag{3.4}
\]

where \(Y,Y'\) are conditionally independent given \(X\).

Here is a direct verification. Write

\[
 X=\sqrt K\,A,\qquad G_s=\sqrt{\alpha K}\,B,\qquad
 G_s'=\sqrt{\alpha K}\,B',
\]

with \(A,B,B'\) independent standard Gaussians, and put

\[
 Z={A/\sqrt\alpha+B\over\sqrt{1+\alpha}},\qquad
 Z'={A/\sqrt\alpha+B'\over\sqrt{1+\alpha}}.
\]

Then

\[
 {W(Y)\over M_s}={e^{-Z^2/2}\over\sqrt{1+\alpha}},\qquad
 \kappa_\lambda(|W(Y)|)
 ={R(Z)\over\sqrt{1+R(Z)^2}},\quad
 R(z)={e^{-z^2/2}\over\sqrt{1+\alpha}}.              \tag{3.5}
\]

On the event

\[
 0\le A\le{\sqrt\alpha\over4},\qquad
 -{1\over4}\le B\le0,\qquad 2\le B'\le{9\over4},     \tag{3.6}
\]

the two values in (3.5) differ by more than \(0.4\). Moreover, using the
minimum of the standard Gaussian density on each displayed interval,

\[
 \mathbb P\{(3.6)\}
 \ge {\sqrt\alpha\over4}\varphi(1/8)\,
      {1\over4}\varphi(1/4)\,
      {1\over4}\varphi(9/4)
 >7\cdot10^{-5}\sqrt\alpha.                          \tag{3.7}
\]

Multiplication by \((1/2)(0.4)^2\) proves (3.4).

Thus even an exact Cheeger minimizer with one perfectly coherent, arbitrarily
long balanced ray has amplitude Dirichlet energy of order
\(\sqrt{s/K}\). Averaging \(\theta\) cannot improve it to the Gibbs-gap
scale \(s/K\). The amplitude obstruction is already present in the
one-dimensional Gaussian equality model; any successful argument must
subtract or otherwise neutralize it rather than bound it by angular
coherence.
