# Angular occupation calculus at the scale \(t=\tau/K\)

## 1. Conclusion

Let ordinary stochastic localization start from a balanced cut, let

\[
 F(t)=\mathbb E I(g_t),\qquad
 b(\tau)={F(\tau/K)\over P_0},\qquad P_0=I(1/2),
\]

and use the \(I(g_t)\)-weighted posterior law

\[
 d\nu_\tau={I(g_{\tau/K})\over F(\tau/K)}d\mathbb P.
\]

Put

\[
 \eta={\sqrt t\,r\over I(g)},\qquad
 \epsilon=1-\eta,\qquad
 m(\tau)=\int\eta^2d\nu_\tau,\qquad
 d(\tau)=1-m(\tau).                                      \tag{1.1}
\]

The symbol \(d\) in this note is the deterministic scalar profile deficit,
not a differential.  If

\[
 \mathfrak a_t={\|P_tD_t\|_{HS}^2\over r_t^2}             \tag{1.2}
\]

is the quadratic-variation density of the active direction, the audited
general angular theorem gives, on

\[
 \mathsf C_{\delta,\alpha}(t)
 =\{g_t\in[\delta,1-\delta],\ \eta_t\ge\alpha\},
\]

\[
 \mathfrak a_t\le {C_{\delta,\alpha}\over t}
     \Omega_\delta(\epsilon_t).                         \tag{1.3}
\]

Here

\[
 \Omega_\delta(s)\lesssim_\delta
 \omega(s):=s^{1/6}\log(e/s)^{1/12},\qquad 0<s\le1,    \tag{1.4}
\]

and \(\omega(0)=0\).  The function \(\omega\) is increasing and concave.
Consequently the factors of \(K\) cancel exactly under \(t=\tau/K\):

\[
\boxed{
 {1\over P_0}\int_{\tau_0/K}^{\tau_1/K}
 \mathbb E\!\left[I(g_t)\mathbf1_{\mathsf C_{\delta,\alpha}(t)}
                    \mathfrak a_t\right]dt
 \le C_{\delta,\alpha}\int_{\tau_0}^{\tau_1}
       b(\tau)\omega(d(\tau)){d\tau\over\tau}.}       \tag{1.5}
\]

This is a dimension-free deterministic weighted occupation lemma.  Its
right side has an exact Jensen reduction stated in Section 4, and is finite
all the way to \(\tau=\infty\) for every positive starting time.  The exact
centroid-defect identity gives a second, independent bound on every stopped
multiplicative window:

\[
\boxed{
 {1\over P_0}\mathbb E\int_{\tau_0/K}^{\sigma}
 I(g_t)\mathfrak a_tdt
 \le {2\over\alpha}\sqrt{\tau_1\over\tau_0}\,
       b(\tau_0)e(\tau_0)
 \le {2\over\alpha}\sqrt{\tau_1\over\tau_0}\,
       b(\tau_0)d(\tau_0),}                          \tag{1.6}
\]

where \(\sigma\le\tau_1/K\) is the first exit from
\(\mathsf C_{\delta,\alpha}\), and

\[
 e(\tau)=\int(1-\eta)d\nu_\tau.                    \tag{1.7}
\]

Thus (1.5) controls arbitrarily long weighted occupation, while (1.6) is
stronger on a fixed multiplicative window whose left endpoint is already
near equality.

There is no converse lower seed in the scalar calculus.  This is sharp in
three separate senses.

1. A Gaussian halfspace is an actual equality model with
   \(\mathfrak a_t=0\), although \(d(\tau)>0\) at every finite time.
2. The heat-flow identity permits every nonnegative scalar remainder to be
   assigned to its nonangular part, with zero angular part; the Gaussian
   model realizes this allocation using one-dimensional observation
   curvature alone.  Hence no nonzero functional of \((b,d)\) can be
   inferred as a lower bound for angular energy from the scalar identities
   alone.
3. The compact perturbations
   \(b_\varepsilon=b_Ge^{\varepsilon\phi}\),
   \(\phi\in C_c^\infty((2,4))\), preserve the initial profile, the final
   perimeter endpoint, and all scalar inequalities while changing every
   non-null local deficit weighting.  The only local weights invisible to
   all such perturbations are endpoint terms.

True T3 extremality may supply a non-scalar phase-rank or focal-incidence
input, but it supplies no scalar-profile lower seed.  No posterior-overlap
claim, KLS bound, or unproved incidence statement is used below.

## 2. Exact scalar profile identities and all endpoints

The ordinary localization identity

\[
 -F'(t)={1\over2t}\mathbb E[I(g_t)\eta_t^2]
\]

becomes

\[
 \boxed{m(\tau)=-{2\tau b'(\tau)\over b(\tau)},
        \qquad d(\tau)=1+{2\tau b'(\tau)\over b(\tau)}.} \tag{2.1}
\]

Define

\[
 h(\tau)=\sqrt\tau\,b(\tau),\qquad
 q(\tau)={b(\tau)d(\tau)\over\sqrt\tau}.            \tag{2.2}
\]

Then, exactly,

\[
 \boxed{h'(\tau)={q(\tau)\over2},\qquad
 {d\over d\tau}\log h(\tau)={d(\tau)\over2\tau}.} \tag{2.3}
\]

For a finite-perimeter balanced cut,

\[
 b(0+)=1,\qquad h(0+)=0,\qquad
 h(\infty)=\lambda={P_\mu(S)\sqrt K\over P_0}.       \tag{2.4}
\]

The two exact finite-endpoint budgets are therefore

\[
 \boxed{
 \int_{\tau_0}^{\tau_1}d(\tau){d\tau\over\tau}
 =2\log{h(\tau_1)\over h(\tau_0)},}                 \tag{2.5}
\]

and

\[
 \boxed{
 h(\tau_1)-h(\tau_0)
 ={1\over2}\int_{\tau_0}^{\tau_1}
       {b(\tau)d(\tau)\over\sqrt\tau}\,d\tau.}   \tag{2.6}
\]

In particular,

\[
 \int_{\tau_0}^{\infty}d(\tau){d\tau\over\tau}
 =2\log{\lambda\over\sqrt{\tau_0}b(\tau_0)},      \tag{2.7}
\]

and

\[
 \lambda-\sqrt{\tau_0}b(\tau_0)
 ={1\over2}\int_{\tau_0}^{\infty}q(\tau)d\tau.    \tag{2.8}
\]

The heat-flow Bernstein identity is

\[
 -q'(\tau)=\mathcal G(\tau)
 +{b(\tau)\over2\tau^{3/2}}
       \int(1-\eta^2)^2d\nu_\tau,                  \tag{2.9}
\]

where \(\mathcal G\ge0\) is the sum of Hessian and observation-curvature
energies.  Jensen gives

\[
 \mathcal E_b(\tau):=-q'(\tau)
 -{b(\tau)d(\tau)^2\over2\tau^{3/2}}\ge0.          \tag{2.10}
\]

In logarithmic time \(x=\log\tau\), write \(d_x=\tau d'(\tau)\).
Direct differentiation gives the useful exact form

\[
 \boxed{\mathcal E_b(\tau)
 ={b(\tau)\over\tau^{3/2}}
       \{d(\tau)(1-d(\tau))-d_x(\tau)\}.}          \tag{2.11}
\]

Thus scalar admissibility is precisely

\[
 0\le d\le1,\qquad d_x\le d(1-d).                  \tag{2.12}
\]

It follows that \(q\) is decreasing.  If \(h(\infty)=\lambda<\infty\),
then (2.7) and (2.12) imply \(d(\tau)\to0\), hence

\[
 m(\tau)\to1,\qquad q(\tau)\to0.                  \tag{2.13}
\]

At the other endpoint, finite second moments give
\(m(\tau)=O(\tau)\), so

\[
 d(0+)=1,\qquad h(\tau)\sim\sqrt\tau,\qquad
 q(\tau)\sim\tau^{-1/2}.                           \tag{2.14}
\]

This singular initial endpoint matters: no angular occupation estimate
derived here can start at \(\tau_0=0\).

The absolute centroid defect remembers a scalar statistic not determined by
\(d\).  Indeed,

\[
 \mathbb E\Delta_{\tau/K}
 =P_0\sqrt{K\over\tau}\,b(\tau)e(\tau),            \tag{2.15}
\]

and, sharply,

\[
 \boxed{1-\sqrt{1-d(\tau)}\le e(\tau)\le d(\tau).} \tag{2.16}
\]

The lower bound is Cauchy--Schwarz, and the upper bound uses
\(\eta^2\le\eta\).  In particular \(d/2\le e\le d\).  This extra scalar
nonuniqueness is one reason that the profile cannot reverse the exact defect
budget.

## 3. The audited pointwise estimate in rescaled variables

At time \(t=\tau/K\), the general angular theorem says

\[
 \|PD\|_{HS}^2\le C_\delta {K^2\over\tau^2}
             \Omega_\delta(\epsilon).               \tag{3.1}
\]

On \(\mathsf C_{\delta,\alpha}\),

\[
 r^2={I(g)^2\eta^2\over t}
 \ge {I_\delta^2\alpha^2K\over\tau},
 \qquad I_\delta=\min_{[\delta,1-\delta]}I>0.
\]

Therefore

\[
 \mathfrak a_{\tau/K}
 \le C_{\delta,\alpha}{K\over\tau}
             \Omega_\delta(\epsilon).              \tag{3.2}
\]

Since \(dt=d\tau/K\), the \(K\)-dependence cancels:

\[
 \mathfrak a_tdt
 \le C_{\delta,\alpha}\Omega_\delta(\epsilon)
        {d\tau\over\tau}.                           \tag{3.3}
\]

The explicit modulus from the audited proof is

\[
 \rho_\delta(s)=C_\delta\left\{\sqrt s+
       (s\sqrt{\log(e/s)})^{1/3}\right\},
 \qquad
 \Omega_\delta(s)=s+\sqrt{\rho_\delta(s)}.         \tag{3.4}
\]

For \(0<s\le1\),

\[
 \Omega_\delta(s)\le C_\delta
       s^{1/6}\log(e/s)^{1/12}=C_\delta\omega(s).   \tag{3.5}
\]

This is the weak \(s^{1/6}\) modulus; no stronger power is silently used.
A direct derivative check shows that \(\omega\) is increasing and concave:
if \(L=\log(e/s)\), then

\[
 {s\omega'(s)\over\omega(s)}={1\over6}-{1\over12L}>0,
\]

and

\[
 {s^2\omega''(s)\over\omega(s)}
 =A^2-A-{1\over12L^2}<0,
 \qquad A={1\over6}-{1\over12L}.                    \tag{3.6}
\]

Finally,

\[
 d(\tau)=\int\epsilon(1+\eta)d\nu_\tau,
 \qquad \int\epsilon d\nu_\tau\le d(\tau).       \tag{3.7}
\]

Concavity and Jensen now give

\[
 \int\Omega_\delta(\epsilon)d\nu_\tau
 \le C_\delta\omega\!\left(\int\epsilon d\nu_\tau\right)
 \le C_\delta\omega(d(\tau)).                      \tag{3.8}
\]

Multiplication by \(I(g)\), expectation, and (3.3) prove (1.5).  Since
\(I(g)\ge I_\delta\) on the central event, its unweighted version is

\[
 \int_{\tau_0/K}^{\tau_1/K}
 \mathbb E[\mathbf1_{\mathsf C_{\delta,\alpha}}\mathfrak a_t]dt
 \le {P_0C_{\delta,\alpha}\over I_\delta}
       \int_{\tau_0}^{\tau_1}b(\tau)\omega(d(\tau)){d\tau\over\tau}.
                                                               \tag{3.9}
\]

## 4. Deterministic integration of the weak modulus

Set

\[
 L_{01}=\int_{\tau_0}^{\tau_1}b(\tau){d\tau\over\tau},
 \qquad
 B_{01}=\int_{\tau_0}^{\tau_1}b(\tau)d(\tau){d\tau\over\tau}. \tag{4.1}
\]

Equation (2.1) gives the exact identity

\[
 \boxed{B_{01}=L_{01}-2\{b(\tau_0)-b(\tau_1)\}.}    \tag{4.2}
\]

In particular \(0\le B_{01}\le L_{01}\).  Applying Jensen to the
probability measure \(b(\tau)d\tau/(\tau L_{01})\) yields

\[
\boxed{
 \int_{\tau_0}^{\tau_1}b(\tau)\omega(d(\tau)){d\tau\over\tau}
 \le L_{01}\omega(B_{01}/L_{01}).}                 \tag{4.3}
\]

Equivalently, when \(B_{01}>0\),

\[
 L_{01}\omega(B_{01}/L_{01})
 =B_{01}^{1/6}L_{01}^{5/6}
       \log(eL_{01}/B_{01})^{1/12}.                 \tag{4.4}
\]

This is the exact deterministic integration of the weak audited modulus.
It does not replace \(1/6\) by a more favorable exponent.

There is also a pure endpoint version.  Since \(h\) is increasing,

\[
 L_{01}\le
 2h(\tau_1)\left({1\over\sqrt{\tau_0}}
                  -{1\over\sqrt{\tau_1}}\right)=:L_*,          \tag{4.5}
\]

and, by (2.5) and monotonicity of \(b\),

\[
 B_{01}\le
 2b(\tau_0)\log{h(\tau_1)\over h(\tau_0)}.         \tag{4.6}
\]

Let

\[
 B_*:=\min\left\{L_*,
 2b(\tau_0)\log{h(\tau_1)\over h(\tau_0)}\right\}. \tag{4.7}
\]

The right side of (4.4) is increasing separately in \(B\) and \(L\) on
\(0<B\le L\).  Hence

\[
 \boxed{
 \int_{\tau_0}^{\tau_1}b\omega(d){d\tau\over\tau}
 \le B_*^{1/6}L_*^{5/6}
       \log(eL_*/B_*)^{1/12}.}                      \tag{4.8}
\]

The right side is interpreted as zero when \(B_*=0\).  Taking
\(\tau_1=\infty\) gives

\[
 L_*={2\lambda\over\sqrt{\tau_0}},\qquad
 B_*=\min\left\{{2\lambda\over\sqrt{\tau_0}},
 2b(\tau_0)\log{\lambda\over\sqrt{\tau_0}b(\tau_0)}\right\}.
                                                               \tag{4.9}
\]

Thus the factor \(b(\tau)\), which is the probability-profile weight of
central posterior states, makes the weak modulus integrable at infinity.
By contrast, the same expression diverges at \(\tau_0=0\), exactly as the
endpoint asymptotics (2.14) predict.

## 5. The exact defect envelope

The exact stochastic identity is

\[
 d\Delta_t=dM_t-\left\{
 {\|P_tD_t\|_{HS}^2\over2r_t}
 +r_tu_t^T(t^{-1}I-A_t)u_t
 +{I(g_t)(1-\eta_t)^2\over2t^{3/2}}
 \right\}dt.                                       \tag{5.1}
\]

Since

\[
 {\|PD\|_{HS}^2\over2r}
 ={I(g)\eta\over2\sqrt t}\,\mathfrak a_t,         \tag{5.2}
\]

the stopped integrated identity gives

\[
 \mathbb E\int_{t_0}^{\sigma}
 {I(g_t)\eta_t\over2\sqrt t}\mathfrak a_tdt
 \le\mathbb E[\mathbf1_{\mathsf C(t_0)}\Delta_{t_0}].          \tag{5.3}
\]

On the stopped window, \(\eta\ge\alpha\) and
\(t\le t_1\).  Therefore

\[
 {1\over P_0}\mathbb E\int_{t_0}^{\sigma}I(g_t)\mathfrak a_tdt
 \le {2\sqrt{t_1}\over\alpha P_0}\mathbb E\Delta_{t_0}.       \tag{5.4}
\]

Substitution of \(t_j=\tau_j/K\) and (2.15) proves (1.6): every factor of
\(K\) cancels.  Combining Sections 3 and 5 gives the useful two-envelope
bound

\[
 \boxed{
 {1\over P_0}\mathbb E\int_{\tau_0/K}^{\sigma}
 I(g_t)\mathfrak a_tdt
 \le C_{\delta,\alpha}\min\left\{
 \int_{\tau_0}^{\tau_1}b\omega(d){d\tau\over\tau},
 \sqrt{\tau_1\over\tau_0}\,b(\tau_0)d(\tau_0)
 \right\}.}                                        \tag{5.5}
\]

The first envelope is useful on a long tail; the second sees a near-equality
left endpoint.  Both are upper bounds.  Neither creates angular energy.

## 6. The heat-flow angular envelope

Let \(s=K/\tau=1/t\), and use the notation of the Bernstein identity.  The
projected-Hessian part of \(\mathcal G\) is

\[
 \mathcal G_{\angle}(\tau)
 ={K^2\over P_0\tau^{7/2}}
   \int\rho_{K/\tau}\|P\nabla_y^2z\|_{HS}^2dy.      \tag{6.1}
\]

The exact posterior dictionary

\[
 P\nabla_y^2z={t^2PD\over I(g)},\qquad
 r={I(g)\eta\over\sqrt t}
\]

gives

\[
 \boxed{
 \mathcal G_{\angle}(\tau)
 ={1\over P_0K\sqrt\tau}\,
   \mathbb E[I(g_{\tau/K})\eta_{\tau/K}^2
             \mathfrak a_{\tau/K}].}               \tag{6.2}
\]

This checks all powers of \(K\) independently of the stochastic
calculation.  Since the angular square is only one summand of
\(\mathcal G\), (2.9)--(2.10) imply

\[
 \boxed{0\le\mathcal G_{\angle}(\tau)
 \le\mathcal G(\tau)\le\mathcal E_b(\tau).}         \tag{6.3}
\]

On the central near-equality event, (3.2), (3.8), and (6.2) also give

\[
 \mathcal G_{\angle}^{\delta,\alpha}(\tau)
 \le C_{\delta,\alpha}{b(\tau)\over\tau^{3/2}}
             \omega(d(\tau)).                       \tag{6.4}
\]

Thus the exact heat-flow two-envelope form is

\[
 \boxed{
 \mathcal G_{\angle}^{\delta,\alpha}(\tau)
 \le\min\left\{\mathcal E_b(\tau),
 C_{\delta,\alpha}{b(\tau)\over\tau^{3/2}}
        \omega(d(\tau))\right\}.}                  \tag{6.5}
\]

Again, both terms are upper budgets.  The first does not distinguish
angular Hessian energy from longitudinal Hessian energy, observation
curvature, or variance of \(\eta^2\).

## 7. Gaussian halfspace test

Take \(\mu=N(0,KI)\), use the centered halfspace, and use the same \(K\) as
the reference scale.  Then

\[
 b_G(\tau)=(1+\tau)^{-1/2},\qquad
 \eta_G^2=m_G={\tau\over1+\tau},\qquad
 d_G={1\over1+\tau}.                                \tag{7.1}
\]

Moreover

\[
 h_G(\tau)=\sqrt{\tau\over1+\tau},\qquad
 \lambda=1,\qquad
 \epsilon_G=1-\sqrt{\tau\over1+\tau}.              \tag{7.2}
\]

The endpoint identities read

\[
 \int_{\tau_0}^{\infty}{d_G(\tau)\over\tau}d\tau
 =\log{1+\tau_0\over\tau_0}
 =2\log{1\over h_G(\tau_0)}.                        \tag{7.3}
\]

The scalar Bernstein quantities are

\[
 q_G(\tau)={1\over\sqrt\tau(1+\tau)^{3/2}},
 \qquad
 \mathcal E_{b_G}(\tau)
 ={2\over\sqrt\tau(1+\tau)^{5/2}}>0.              \tag{7.4}
\]

But the active direction is constant and

\[
 PD=0,\qquad \mathfrak a_t=0,\qquad
 \mathcal G_{\angle}=0.                             \tag{7.5}
\]

All of (7.4) is one-dimensional observation curvature.  In particular,

\[
 \int_{\tau_0}^{\infty}b_G(\tau)\omega(d_G(\tau)){d\tau\over\tau}>0
\]

for every \(\tau_0>0\), while the angular occupation is exactly zero.  At
infinity its integrand is
\(O(\tau^{-5/3}\log(\tau)^{1/12})\), so it is integrable.  At zero it is
asymptotic to \(1/\tau\), so it diverges.  This simultaneously checks the
tail conclusion and the excluded endpoint in Section 4.

More generally, for a Gaussian reference whose tail coefficient is
\(\lambda>0\),

\[
 b_{G,\lambda}(\tau)={\lambda\over\sqrt{\lambda^2+\tau}},\qquad
 d_{G,\lambda}(\tau)={\lambda^2\over\lambda^2+\tau},             \tag{7.6}
\]

and

\[
 q_{G,\lambda}(\tau)
 ={\lambda^3\over\sqrt\tau(\lambda^2+\tau)^{3/2}},
 \qquad
 \mathcal E_{G,\lambda}(\tau)
 ={2\lambda^3\over\sqrt\tau(\lambda^2+\tau)^{5/2}}.           \tag{7.7}
\]

The angular energy is still zero.  Thus changing the perimeter endpoint
does not by itself create a scalar angular seed.

## 8. Compact scalar perturbations and the sharp lower-seed no-go

Let \(\phi\in C_c^\infty((2,4))\) be nonzero and put

\[
 b_\varepsilon(\tau)=b_G(\tau)e^{\varepsilon\phi(\tau)}.        \tag{8.1}
\]

For all sufficiently small positive and negative \(\varepsilon\),
\(b_\varepsilon\) is scalar-admissible.  Indeed,

\[
 m_\varepsilon=m_G-2\varepsilon\tau\phi'(\tau),
 \qquad
 d_\varepsilon=d_G+2\varepsilon\tau\phi'(\tau),                \tag{8.2}
\]

and on \([2,4]\) the Gaussian inequalities

\[
 b_G'<0,\qquad 0<d_G<1,\qquad
 d_{G,x}<d_G(1-d_G)                                  \tag{8.3}
\]

have strict uniform margins.  Outside \([2,4]\), the perturbed and
Gaussian profiles agree.  Therefore

\[
 b_\varepsilon(0+)=1,\qquad
 \sqrt\tau b_\varepsilon(\tau)\longrightarrow1,                 \tag{8.4}
\]

and every seed value below time \(2\) is unchanged.  Also

\[
 \int_0^\infty(d_\varepsilon-d_G){d\tau\over\tau}
 =2\varepsilon\int_0^\infty\phi'(\tau)d\tau=0,                  \tag{8.5}
\]

where the first equality is understood on any interval containing the
support.  Thus the exact log-endpoint budget is preserved.

There is an algebraic zero-angular allocation for every such scalar
profile at the level of all scalar profile and heat-flow identities used
in the lower-seed question.  Set
\(\eta^2=m_\varepsilon\) deterministically under the formal boundary law,
so that

\[
 \operatorname {Var}_{\nu_\tau}(\eta^2)=0,
\]

set \(\mathcal G_{\angle}=0\), and assign

\[
 \mathcal G_{\rm nonang}(\tau)=\mathcal E_{b_\varepsilon}(\tau)\ge0.
                                                               \tag{8.6}
\]

Then (2.9) holds exactly, because

\[
 -q'_{b_\varepsilon}
 =\mathcal E_{b_\varepsilon}
 +{b_\varepsilon d_\varepsilon^2\over2\tau^{3/2}}.             \tag{8.7}
\]

This is a formal scalar counterprofile, not a claim that every
\(b_\varepsilon\) is realized by a log-concave measure or by a pathwise
stochastic localization.  Its force is precise: any proof using only
\(b,d\), their endpoints, profile monotonicity, and the scalar
sum-of-squares identity is compatible with zero angular energy.  Therefore,
if \(\mathcal W[b,d]\ge0\) is any
deterministic scalar functional and those scalar facts alone implied

\[
 \int\mathcal G_{\angle}\ge\mathcal W[b,d],          \tag{8.8}
\]

then \(\mathcal W\) would have to vanish on every scalar-admissible
counterprofile.  There is no nontrivial scalar-only lower seed.

The exact stochastic defect identity does not reverse this conclusion.
After expectation it controls only the sum of angular dissipation, active
covariance defect, and the scalar \((1-\eta)^2\) square.  It gives no lower
bound on the first summand, and (2.16) shows that \(d\) does not even
determine its total initial budget \(\mathbb E\Delta\).  The Gaussian
halfspace is an actual pathwise model in which the angular summand is zero
and the other nonnegative summands pay the defect drop.

The compact perturbations also characterize the exceptional local
weightings exactly.  Put \(x=\log\tau\), \(u=\log b\).  Then

\[
 d=1+2u_x.                                           \tag{8.9}
\]

For a local weighted functional

\[
 \mathcal J_w[b]=\int b(\tau)w(d(\tau)){d\tau\over\tau},        \tag{8.10}
\]

the perturbation \(u\mapsto u+\varepsilon\psi\),
\(\psi\in C_c^\infty\), gives

\[
 {d\over d\varepsilon}\mathcal J_w[b e^{\varepsilon\psi}]
 \bigg|_{\varepsilon=0}
 =\int b\{w(d)\psi+2w'(d)\psi_x\}\,dx.             \tag{8.11}
\]

If this variation vanishes for every compact perturbation of every strict
admissible profile, its Euler equation is

\[
 w(d)+(1-d)w'(d)-2w''(d)d_x=0.                       \tag{8.12}
\]

Since \(d_x\) varies in an open interval under (2.12), one must have
\(w''=0\), and then (8.12) forces

\[
 w(d)=c(1-d).                                        \tag{8.13}
\]

But this is only an endpoint term:

\[
 \int b(1-d){d\tau\over\tau}
 =\int bm\,dx=2\{b(\tau_0)-b(\tau_1)\}.            \tag{8.14}
\]

Every other \(C^2\) local deficit weighting, including
\(w(d)=\omega(d)\), changes under some compactly supported perturbation
with the same endpoints.  Without the factor \(b\), the analogous null
weights are affine functions of \(d\); their \(d\)-part is exactly the
log-endpoint term (2.5).  These are the only scalar null Lagrangians, and
none measures angular energy.

## 9. Consequence for a T3 extremizer

T3 extremality supplies genuine geometric information: signed-distance
rays, constant-mass separated tails, and, in the bad-scale reduction,
high-rank phase information.  None of those facts is a scalar consequence
of \((b,d)\).  Conversely, the normalized scalar endpoints record only

\[
 b(0+)=1,\qquad
 \sqrt\tau b(\tau)\to {P_\mu(S)\sqrt K\over P_0}.    \tag{9.1}
\]

The Gaussian linear/halfspace branch is an exact equality model for the
posterior centroid and angular identities and has zero angular energy.
Compact perturbations preserve (9.1), preserve any finite initial seed
interval, and preserve the exact scalar inequalities.  Thus the presently
proved consequences of T3 extremality do not turn any weighting of the
scalar profile deficit into a positive angular lower bound.

A successful lower seed must use an additional non-scalar statement, for
example a proved inequality of the schematic form

\[
 \text{phase rank or focal incidence}
 \ \Longrightarrow\ 
 \int\mathcal G_{\angle}>0.                          \tag{9.2}
\]

The current Fisher-rank and transport-ray reductions suggest such an
input, but they do not prove it: converting phase diversity into common
posterior mass or interface capacity is exactly the unresolved overlap and
incidence problem.  Importing that conversion would be circular at KLS
strength.  The rigorous output of the present calculus is therefore the
upper occupation bounds (1.5), (1.6), and (6.5), together with the sharp
scalar-only lower-seed no-go of Section 8.
