# Chapter Plan — *Gaussian Natural Gradient Flow for Variational Inference* (FoCM)

Planning document. Nothing in `manuscript/` is edited until this plan is confirmed.
Companion to (and superseding, for structure) Codex's `manuscript-plan.md`.

---

## 0. Locked decisions (INSIGHT collection)

Recorded in the author's own words / confirmed choices:

- **[INSIGHT: thesis]** "I want the readers to remember the **three-stage phenomena** of Gaussian
  approximate Fisher–Rao gradient flow (both continuous-time dynamics and deterministic or
  stochastic discretization schemes); the general affine invariant metric is just an extension
  point."
- **[INSIGHT: framing]** VI-method paper: *"Gaussian natural gradient flow for variational
  inference"*. Title follows `natural-gradient.tex` ("Gaussian Natural Gradient Flows for
  Variational Inference"), not the rates-first Codex title.
- **[INSIGHT: style]** Writing follows `natural-gradient.tex` — the author's own hand-written
  manuscript — not the audited-note voice of the improved notes and not Codex's draft voice.
- **[INSIGHT: S4↔S3]** Section 4 "can refer to the exact theorem/lemma or other stuff in
  Section 3, but should not just call the names" — i.e., restate the content it uses
  (contraction factors, radii, stepsize windows) in place, with a `\cref` pointer.
- **[INSIGHT: structure]** Five main sections as specified by the author; geometry
  classification is the closing Section 5; non-log-concave material lives inside Section 2 with
  a *brief* instability example motivating clipping; BW warm start / quadratic rescue appear at
  the end of the log-concave part of Section 2 as corollary-level "burn-in removal"; the
  BW comparison is at most a remark.
- **[INSIGHT: venue]** Foundations of Computational Mathematics. Submission PDF first; EMS
  template only after acceptance — keep the portable `article` setup.
- Sections 1–3 are entirely deterministic; Section 4 is the stochastic theory; WFR is
  explicitly a second paper (cut).
- **[Confirmed 2026-07-22]** Small focused numerics section (§6, ~3–4 composite figures tied
  to theorems; protocols in an appendix); the legacy *stochastic* clipped-KL theorem is
  **dropped** (mentioned only as an open direction in the Discussion); notation follows
  `natural-gradient.tex` ($N_\theta$, $\rho_\post$, $\mathcal E$, $\Delta t$); drafting is a
  **fresh rewrite** of `sections/*.tex` in the author's voice with the Codex files kept aside
  as quarry.

---

## 1. Style contract

### 1.1 Voice (from `natural-gradient.tex`)

- First-person plural, walking narrative. Nearly every result gets exactly **one lead-in
  sentence**: "We proceed to show the convergence theory of \eqref{upd:Riem} under
  Assumption \ref{assump:logconcave-smooth}." / "We first show the spectral bounds of the
  covariance along \eqref{ODEs:NG}."
- Motivation is **one pragmatic sentence**, not a paragraph: "…is geometrically natural but
  computationally expensive… for the ease in computation, we develop the following
  discretization scheme with KL divergence."
- Hypotheses inline in statements ("Under Assumption \ref{...}, if the initial covariance
  satisfies …, then …"). Proofs are direct computations with short connective phrases
  ("This gives", "by Gronwall's inequality"), no proof-sketch meta-commentary.
- No bold run-in "contribution paragraph" headers in the body; no "Purpose of the note"-style
  meta-discourse; no defensive hedging blocks. Scope caveats become short remarks.

### 1.2 Architecture lessons adopted from the two model papers

From Lambert–Chewi–Bach–Bonnabel–Rigollet (2205.15902) and Diao–Balasubramanian–Chewi–Salim
(2304.05398):

1. **Display the central objects in the introduction**: the KL objective and the boxed NG flow
   ODEs appear in §1 before any formal setup, plus the boxed three-stage decomposition
   `covariance burn-in + global localization + local convergence` as the paper's organizing
   device.
2. **Contributions as flowing prose with numbered forward-pointers** ("we show in
   Theorem 3.2 that …"), no informal-theorem environments, no bullet lists, no rates in the
   intro prose — quantitative statements wait for the numbered theorems.
3. **Economy of numbered main-text statements.** The notes carry ~200 formal items; the paper's
   main text should number roughly 40–50, with headline theorems clearly identifiable
   (Lambert numbers 5 statements in 12 pages; we are denser but must stay selective).
   Technical/supporting lemmas go to appendices.
4. **Caveats live in short post-theorem remarks**, never inside theorem statements and never in
   a limitations dump ("The upper bound … is notationally convenient for our proof but not
   necessary" is the register).
5. **Counterexamples: one honest sentence in main text + a statement-level
   corollary, full construction in a dedicated appendix.**
6. **Complexity embedded in the theorem** Diao-style: clean rate first, then "In particular,
   with $\Delta t \asymp …$ and $N \gtrsim …$, we obtain $\varepsilon$" in the same statement
   (or an immediately following corollary), not a separate scattered corollary.
7. **One pseudocode box** for the stochastic pipeline in implementable $(m, C)$ coordinates.
8. **Related work**: thematic paragraphs at the end of §1, each ending with a positioning
   sentence; the Wasserstein/BW line is treated respectfully as "a different geometry"
   (mirroring how Lambert et al. treat Fisher–Rao), engaged quantitatively only where we use
   it (warm start, nonconvex comparison remark).
9. **Conclusion short and forward-leaning**: recap in two paragraphs, concrete open problems,
   no re-litigation of limitations.

### 1.3 Notation (confirmed)

Keep `natural-gradient.tex` conventions throughout: $N_\theta$ (dimension),
$\rho_\post$, $\mathcal E$, $\mathcal A$, $a=(m,C)$, $\Delta t$, $\PG$. The improved notes'
whitened quantities are imported with their names ($\alpha_\star, \beta_\star, \kappa_\star,
\lambda_{0,\star}, \Gamma$, $\Delta(a)$). One notation paragraph at the end of §1; no notation
table.

---

## 2. Section-by-section plan

Page estimates are for the main text; FoCM has no hard limit but ~50 pp main text +
appendices is the target envelope.

### §1 Introduction (~5 pp)

**Purpose.** Frame the paper as *the* convergence theory of Gaussian natural gradient VI;
fix the reader's takeaway (three-stage phenomenon) before any notation.

**Reader takeaway.** Fisher–Rao Gaussian VI converges in three stages — covariance burn-in,
global localization at the intrinsic condition number, local convergence at a score-operator
spectral scale — in continuous time and under both deterministic and stochastic
discretizations, with matching lower bounds at each stage.

| unit | content | source |
|---|---|---|
| opening | sampling → VI → Gaussian family → KL objective (displayed) | `natural-gradient.tex` §1, furnished |
| flow | **boxed** NG flow ODEs `ODEs:NG`; FR = natural gradient; affine invariance in one paragraph (seeds §5) | `natural-gradient.tex` + det-note whitening vocabulary |
| three-stage box | boxed `burn-in + global localization + local convergence` display + one paragraph per stage, in plain language | det note l.231-239; stoch note intro |
| contributions | prose paragraphs with forward pointers: (i) continuous global rate + spiral sharpness; (ii) two discretizations, global rates, $\kappa/\Delta t$ sharpness, burn-in removal; (iii) non-log-concave stationarity + clipped scheme; (iv) local spectral sandwich + Gaussian-core region + convex-ridge sharpness; (v) deterministic three-stage complexity; (vi) stochastic STL/Price–Hessian theory: floors, decreasing-step, stochastic three-stage; (vii) affine-invariant classification as extension | new prose; raw material in Codex `01-introduction.tex` (three-bottleneck framing) — rewrite in author's voice |
| related work | Wasserstein/BW VI (Lambert, Diao); natural-gradient & information geometry (Rao, Amari, Chen et al.); STL (Roeder et al.); nonconjugate/constrained NGVI (Sun et al.); positioning sentence each | `01-introduction.tex` §1.3 + `natural-gradient.tex` citations |
| organization + notation | one short paragraph each | new |

Cut here: Codex's abstract framing ("optimizer whitening separates three rate-limiting
mechanisms" is kept as an idea, rewritten); the intro carries **no** theorem environments.

### §2 Global convergence and discretization (~14 pp)

**Purpose.** The full global deterministic theory: continuous flow, two discretization
schemes, sharpness at both levels, burn-in removal, and the non-log-concave boundary.

**Reader takeaway.** Global convergence costs
$\logp(1/(\beta\lambda_{0,\min})) + \kappa\text{-scale localization}$, this is sharp in both
continuous and fixed-step senses, the burn-in term is removable by transport or one Hessian
query, and beyond log-concavity the unclipped schemes genuinely fail while a clipped KL step
retains an $O(1/N)$ stationarity guarantee.

**2.1 Setup and the flow.** Assumption `assump:logconcave-smooth`; FR metric on Gaussians;
Riemannian gradient and flow identification; exact dissipation
$\frac{d}{dt}\Delta(a_t) = -\|\grad\mathcal E\|_{a_t}^2$; the two basic inequalities
(Wasserstein PL on Gaussians + FR/W2 comparison).
*Sources:* `natural-gradient.tex` §2 opening (voice), det note `lem:metric-dissipation`,
`lem:W2-PL`+`lem:FR-W2` (state once; resolves the red TODO tcolorbox — the corrected
constants come from the det note).

**2.2 Continuous-time global convergence.** Precision-Duhamel + covariance bootstrap
(`lem:precision-duhamel`, `lem:cont-cov`); optimizer whitening and intrinsic conditioning
($\alpha_\star,\beta_\star,\kappa_\star$, affine invariance of $\kappa_\star$,
$\beta_\star\lambda_{0,\star}$); **Theorem (headline):** affine-invariant global rate
$\Delta(a_t)\le\Delta_0[1+\beta_\star\lambda_{0,\star}(e^t-1)]^{-1/\kappa_\star}$ with
embedded hitting time $\logp(1/(\beta_\star\lambda_{0,\star})) + \kappa_\star\logp(\Delta_0/\delta)$.
Post-theorem remarks: `rem:kappa-star-scope`; the 1D exact example `rem:FR-init-cov`
(initial-covariance dependence is not a proof artifact — final sentence rewritten to point at
the BW warm start in §2.6 instead of the cut WFR section).
Sharpness: **Corollary** ($N_\theta\ge2$: localization time $\Theta(\kappa_\star)$ under the
spectral initialization class) + one honest sentence pointing to the spiral construction in
Appendix B.
*Sources:* det note l.700–1715 (`thm:cont-global` whitened form supersedes the stoch note's
original-coordinate version; `thm:spiral-sharp`, `cor:cont-sharp`).

**2.3 Discretization with Riemannian distance.** Forward–backward proximal derivation of
`upd:Riem`; equivalence with the retraction step (`prop:upd-map-Riem`); one-sentence remark
that the covariance block is the SPD exponential but the joint map is a retraction (Codex
polish); covariance bootstrap `lem:Riem-cov` (floor $1/(2\beta)$ after
$O(\frac1{\Delta t}\logp\frac1{\beta\lambda_0})$ steps); one-step descent `lem:Riem-descent`;
**Theorem:** global hitting time
$N = O(\frac1{\Delta t}\logp\frac1{\beta\lambda_0} + \frac{\kappa}{\Delta t}\logp\frac{\Delta_0}{\delta})$
for $\Delta t \le 1/(2\beta\lambda_{\max})$.
*Sources:* derivation from `natural-gradient.tex` §2.1 (keep verbatim modulo notation);
rates from det note `thm:Riem-global` (supersedes `thm:conv-Riem`).

**2.4 Discretization with KL divergence.** Bregman/mirror identification of `upd:KL`
(load-bearing, keep in main text); covariance bootstrap `lem:KL-cov` with the *exact* solvable
envelope, valid for **every** $\Delta t>0$ — contrast sentence vs the Riemannian stepsize
restriction; forward-KL one-step descent `prop:KL-onestep`; **Theorem:** same-form global
hitting time. Comparison remark: matrix exponential vs rational resolvent; stability ranges.
*Sources:* derivation from `natural-gradient.tex` §2.2; rates from det note
`thm:KL-global` (supersedes `thm:conv-KL`).

**2.5 Sharpness of fixed-step localization.** **Theorem statement in main text:** for
$\Delta t=\gamma/\kappa$, both schemes need $\Omega(\kappa/\Delta t)$ iterations for
constant-factor reduction, even after burn-in, with dimension-free Hessian-Lipschitz constant;
+ corollary ($\kappa/\Delta t$ factor in 2.3/2.4 is sharp; $\Theta(\kappa^2)$ at
$\Delta t=\Theta(1/\kappa)$). Mechanism remark (`rem:why-slow`) + scope remark (fixed step
only; adaptive/implicit schemes open). Bump train → Appendix C.
*Sources:* det note `thm:disc-global-sharp`, `cor:disc-global-sharp`.

**2.6 Removing the covariance burn-in.** Two corollary-level devices, presented after the
rates they improve: (a) **BW transport warm start** — additive vs multiplicative covariance
growth (one contrast sentence), `lem:W-cov-bootstrap`, continuous and discrete switch
theorems (`thm:cont-W-to-FR`, `thm:disc-W-to-FR`); (b) **quadratic rescue** — one
gradient–Hessian query lands $C_0=H_{\rm in}^{-1}$ on the curvature band
(`def:rescue`, `lem:rescue-band`), exact recovery of Gaussian targets in one query
(`thm:gauss-recovery`), burn-in-free complexity (`thm:disc-rescue-to-FR`). Oracle-tradeoff
remark (`rem:two-initializers`) + caveat remark: these remove only the burn-in stage, not the
$\kappa/\Delta t$ localization factor nor the local stiffness.
*Sources:* det note §Transport bootstrap + stoch note §8 (rescue, sharper).

**2.7 Non-log-concave targets.** In order: (a) continuous-time stationarity
$\min_{t\le T}\|\grad\mathcal E\|^2 \le (\mathcal E(a_0)-\mathcal E(a_\star))/T$ (kept from
`natural-gradient.tex`, trivial proof in place); (b) *brief* instability: the smooth double
well $V_R(x)=x^2/2-2R^2\log\cosh(x/R)$ with $|V_R''|\le1$ — condensed statement(s) of the
Riemannian covariance cascade and the KL resolvent pole (one proposition each, short proofs
possibly merged or deferred to Appendix G), one sentence of moral: a global Hessian bound
cannot replace covariance control; (c) **clipped KL/Bregman scheme**: spectral clipping =
log-det Bregman projection, split stationarity measure $\mathsf S_n$, **Theorem:** descent
$\mathcal E(a_{n+1})\le\mathcal E(a_n)-\frac1h\mathsf S_n$ and
$\min_n \mathsf S_n \le \frac hN(\mathcal E(a_0)-\mathcal E(a_N))$ for
$h\le 1/(2\beta\lambda_+)$; **Corollary:** full successive-Gaussian KL with the
$\lambda_+/\lambda_-$ factor; scope remark (constraint-pinned iterates). BW nonconvex
comparison (Diao Thm 5.5) compressed to a **remark**. Full proofs → Appendix G.
*Sources:* Codex `06-nonconvex.tex` + `appendices/D` (its clipped theorem is the *repaired*
version of `natural-gradient.tex`'s and supersedes it; the ng proof's unhandled clipping
interaction is fixed by the Bregman-projection argument).

### §3 Local convergence rates (~11 pp)

**Purpose.** The local theory for flow + both schemes via the Gaussian core, its sharpness,
and the deterministic three-stage headline theorem.

**Reader takeaway.** Near the optimizer the rate is governed by a self-adjoint score operator
with spectral gap in $[\max\{\alpha_\star,\frac1{4+\Gamma}\},\ \min\{\beta_\star,\frac{4+\Gamma}2\}]$,
$\Gamma=O(1+\sqrt{N_\theta}\log\kappa_\star)$; the certified local region is *exact* for
Gaussian targets; the dimension dependence of the gap is real (convex ridge); and the full
picture assembles into the three-stage complexity.

**3.1 Whitened coordinates and the linearized generator.** Whitened setup; first/second-
Hermite score operators $T, T^*, H$; diagonal-mode bounds and
$\Gamma=\min\{\beta_\star-1,\ \sqrt{N_\theta}(4\sqrt{2t_0}+4t_0+4)\}=O(1+\sqrt{N_\theta}\log\kappa_\star)$;
**Proposition (headline):** spectral sandwich
$\max\{\alpha_\star,\frac1{4+\Gamma}\}\le\gamma_\star\le\Lambda_\star\le\min\{\beta_\star,\frac{4+\Gamma}2\}$,
equality picture for Gaussian targets; linearized per-step factors for both maps
(Riemannian $1-h/(4+\Gamma)$; KL $1-h/(4+2h+\Gamma)$ in the $\Delta t$-weighted norm — keep
the weighted-norm self-adjointization trick visible, it is the author's device).
*Sources:* det note l.3034–3370 / stoch note §7 (identical); operator framework skeleton from
`natural-gradient.tex` §Local (superseded rates dropped).

**3.2 The Gaussian core and the certified region.** Entropy-enhanced coercivity + covariance
band; energy sublevel $U_\delta$, hull, non-Gaussian moduli; exact Gaussian-core dissipation
`lem:gauss-core`; non-Gaussian modulus $K_{ng}(\delta)$ (statement only; derivation →
Appendix A); **Theorem:** continuous local convergence at rate
$\rho_{\rm core}(\delta)$, with corollaries: universal level $\to 2/(4+\Gamma)$, and
Gaussian targets attain the *exact* threshold (dimension-free, sharpness benchmark
$\Delta^\#_G(\rho)=\phi(\rho/2)$).
*Sources:* det note §Local (Gaussian core) = stoch note §7.

**3.3 Discrete local convergence.** Exact Gaussian maps: KL scalar map + closed-form
invariant interval (`prop:kl-gauss`, `prop:kl-sharp`); Riemannian exact map + no-overshoot
interval (`prop:riem-gauss`); **Theorem:** general-target local contraction for both schemes
at $(1-\Delta t\,\gamma_\bullet/2)^n$ under the constructive radius condition, with the
honest remark that $K_{\rm disc}$ is constructive-not-closed-form; pointwise one-step
contraction lemma (`lem:local-map-contraction`, stated here because §4 reuses it — per
[INSIGHT: S4↔S3] its content is restated there).
*Sources:* det note §Local discrete + stoch note `lem:local-map-contraction`.

**3.4 Online entry criteria.** Optimizer-free BW-residual gate
($R_{BW}(a)^2\le2\alpha\delta_{\rm loc}$ certifies entry) + exact Gaussian spectral gates;
one remark on stochastic caveat (forward pointer to §4).
*Sources:* det note §entry (`prop:entry-residual`, `prop:entry-gauss`).

**3.5 Sharpness: no dimension-free $\log\kappa_\star$ rate.** **Theorem statement in main
text (integration-ready two-sided form):** a smooth convex-ridge family with
$N_\theta=\Theta(\kappa_\star^2)$, uniformly Hessian-Lipschitz ($L_H=64$), satisfying exact
whitening, has $\tau_H^2=\Theta(\kappa_\star^{1/2})$ and
$\gamma_\star=\Theta(\kappa_\star^{-1/2})$ — saturating the $\beta_\star-1$ branch of
$\Gamma$ and ruling out every dimension-free logarithmic local rate. Scope remark: does not
settle fixed dimension or sharpness of the $\sqrt{N_\theta}\log\kappa_\star$ branch.
Construction → Appendix D.
*Sources:* `local-log-kappa-counterexample.tex` `thm:main` + the commented-out
integration-ready corollary (l.1661–1670) + `rem-scaling-scope`.

**3.6 Deterministic three-stage complexity (paper headline).** **Theorem:** the boxed
three-stage bound for the flow and both schemes,
$T = O(\logp\frac1{\beta_\star\lambda_{0,\star}} + \kappa_\star\logp\frac{\Delta_0}{\delta_c} + (4+\Gamma)\log\frac1\varepsilon)$
and the $/\Delta t$ discrete versions; **Corollary:** burn-in-free version under rescue or
transport warm start; restart remark. Each stage's sharpness cross-referenced in one
sentence each (spiral; bump train; convex ridge).
*Sources:* det note `thm:det-three-stage` + stoch note `cor:burnin-free`, `rem:local-restart`.

### §4 Stochastic algorithms (~11 pp)

**Purpose.** The STL/Price–Hessian stochastic theory reproducing the three-stage picture in
expectation, with explicit floors, their removal, and the stochastic three-stage corollary.

**Reader takeaway.** With one deterministic quadratic rescue and minibatch Price/Hessian STL
iterations, both schemes contract globally in expectation to explicit variance floors
(two assumption regimes), contract locally at the deterministic rate in the killed/non-exit
sense with a floor $\propto L_{H,\star}^2/B$, admit floor-free decreasing-step convergence,
and are exactly calibrated on Gaussian targets — the same three stages, now with floors.

**4.1 Oracle, estimators, and the algorithm.** Price/Hessian oracle; why STL: the plain
single-sample estimators' variance asymmetry (gradient noise $\propto\beta^2\lambda_{\max}N_\theta$
vs a.s.-bounded Hessian noise — content of `natural-gradient.tex` `lem:var-bound-stoch`,
compressed to motivation; the plain-estimator theorems themselves are cut, remark
`rem:price` notes the gradient-only estimator is out of scope); STL mean estimator + sampled
Hessian with a.s. spectral band; quadratic rescue in its stochastic role + rescue-energy
lemma; **pathwise covariance bands with no clipping/stopping** (`lem:stoch-band`; Riemannian
band for $\Delta t\le1/\kappa$, KL band for every $\Delta t$); **Algorithm 1** pseudocode box
in $(m,C)$ coordinates. Fix the `natural-gradient.tex` STL $\hat K_n$ typo (l.829) when
importing the estimator display.
*Sources:* stoch note §10 (`def:oracle`, `alg:unified`, `lem:rescue-energy`,
`lem:stoch-band`); motivation from `natural-gradient.tex`.

**4.2 Intrinsic STL variance.** $\Psi(a)=\mathbb E\|C^{1/2}(\nabla^2V(X)-A)C^{1/2}\|_F^2$;
noise bound $\le 2G^2+\tfrac32\Psi$ (÷$B$ for minibatch); **two regimes**: (a)
assumption-minimal self-bounding $\Psi\le 2\beta\Lambda N_\theta+2\beta\Lambda G^2$; (b)
Hessian-Lipschitz $\Psi\le L_H^2\,{\rm Tr}(C^2){\rm Tr}(C)\le N_\theta^2\bar L_H^2$
(anisotropic form stated before the worst case); $\Psi\equiv0$ for Gaussian targets.
*Sources:* stoch note §11.

**4.3 Global convergence with floors.** One-step lemmas for both schemes (restating the §2
descent structure they perturb, per [INSIGHT: S4↔S3]); **Theorem pair** (or 2×2 presented as
two theorems + a comparison display): Riemannian and KL, Hessian-Lipschitz regime
($\widetilde O(\kappa^2)$ to floor $O(\kappa N_\theta^2\bar L_H^2/B)$) and assumption-minimal
regime ($\widetilde O(\kappa^3)$ single-sample; $\widetilde O(\kappa^2)$ rounds at batch
$\Theta(\kappa)$; floor $O(\kappa N_\theta)$); batch/stepsize-window remarks;
`lem:KL-logtangent` (KL resolvent *is* an exponential retraction — the structural bridge
that lets both schemes share the perturbation analysis; keep visible).
*Sources:* stoch note §12–13 (`thm:Riem-Hlip`, `thm:Riem-selfbound`, `thm:KL-Hlip`,
`thm:KL-selfbound`).

**4.4 Local convergence.** Whitened local modulus $L_{H,\star}$; stochastically admissible
level; local variance decomposition ($V_0=3N_\theta^2L_{H,\star}^2$,
$V_1=24+6N_\theta L_{H,\star}^2$); one-step perturbation of the local maps (restate the §3.3
pointwise contraction content in place); **Theorems:** killed/non-exit mean-square
contraction for both schemes at per-step factor $1-\frac{\Delta t}{2(4+\Gamma)}$
(resp. KL weighted variant) with floor $O(\Delta t(4+\Gamma)N_\theta^2L_{H,\star}^2/B)$;
**Proposition:** finite-horizon exit probability, deep-entry + horizon-batch condition;
honest-scope remark (killed not stopped; STL noise has unbounded support, so infinite-horizon
containment is not claimed); **headline complexity remark:**
$N_{\rm loc}=O((4+\Gamma)^2\max\{1,V_1/B\}\log(1/\varepsilon))$, i.e.
$O((1+\sqrt{N_\theta}\log\kappa)^2\log\frac1\varepsilon)$ under the universal certificate,
$\kappa$-independent under uniform local conditioning; no intermediate $O(\log\kappa)$
theorem claimed. **Corollary:** stochastic global-to-local three-stage pipeline
(`cor:stoch-three-stage`) — the stochastic mirror of §3.6.
*Sources:* stoch note §14 (`thm:stoch-local-Riem`, `thm:stoch-local-KL`,
`prop:finite-horizon`, `rem:local-condnum`, `rem:stoch-local-scope`, `cor:stoch-three-stage`).

**4.5 Removing the floors.** (a) **Decreasing stepsizes:** floor-free $O(1/N)$ under
Hessian-Lipschitzness ($\Delta t_n=8\kappa/(n+n_0)$), $K=0$ for Gaussian targets; (b)
**Gaussian calibration:** rescue recovers the optimizer in one query; matched covariance
$C_0=C_\star$ gives pathwise-zero noise and deterministic contraction; generic $C_0$ keeps
mean noise until the covariance matches (`cor:stoch-local-gauss`, `prop:gauss-nofloor`).
*Sources:* stoch note §15–16.

**4.6 Complexity summary.** One table: regimes × schemes × (contraction, floor, iteration
scale, oracle pairs), constant-step vs decreasing-step; mirrors the boxed three-stage
display one last time.
*Sources:* stoch note §17 (`tab:complexity`).

### §5 General affine-invariant geometry (~5 pp)

**Purpose.** The extension point: classify all affine-invariant metrics on Gaussian space and
show Fisher–Rao is the unique balanced one — explaining, post hoc, the choice made in §2–4.

**Reader takeaway.** Up to scale, affine-invariant Gaussian metrics form the three-parameter
family $(\eta,\omega,\tau)$; on a Gaussian target the parameters precondition the mean, trace,
and traceless covariance modes at exact rates $1/\eta$, $1/(2(\omega+N_\theta\tau))$,
$1/(2\omega)$; Fisher–Rao $(1,\tfrac12,0)$ is the unique metric balancing all three.

| unit | content | source |
|---|---|---|
| 5.1 classification | **Theorem** ($N_\theta\ge2$): the $(\eta,\omega,\tau)$ family with positivity constraints, uniqueness; $N_\theta=1$ identifiability remark; Schur-lemma proof → Appendix H | `natural-gradient.tex` `thm:affine-invariant-metrics` (author's derivation voice) upgraded by Codex `thm:geometry-classification` (adds $d\ge2$, uniqueness; fixes the ng proof's SPD-classification typo at l.2826) |
| 5.2 general flow + modal rates | general $(\eta,\omega,\tau)$ KL gradient-flow ODEs; **exact** trace–traceless decomposition (not merely linearized); **Theorem:** exact Gaussian modal rates $\gamma_{\rm mean}=1/\eta$, $\gamma_{\rm traceless}=1/(2\omega)$, $\gamma_{\rm trace}=1/(2(\omega+N_\theta\tau))$ | `natural-gradient.tex` `ODEs:affine-invariant` + Codex `thm:geometry-general-flow`, `thm:geometry-gaussian-modes` |
| 5.3 balance | **Theorem:** rates equal iff $(\eta,\omega,\tau)=c(1,\tfrac12,0)$ — Fisher–Rao uniquely balanced; scope remark (exact linearization rates at the Gaussian optimizer, *not* general-target global rates; an unbalanced metric can win on a volume-dominated task) | Codex `thm:geometry-balanced`, `rem:geometry-balance-scope` |
| 5.4 general discretizations | the $(\eta,\omega,\tau)$ Riemannian and KL-type schemes, recovering `upd:Riem`/`upd:KL` at $(1,\tfrac12,0)$; keep the author's "not the exact update of a proximal algorithm" caveat | `natural-gradient.tex` `upd:affine-invariant-Riem` + KL variant |

### §6 Numerical illustrations (~3 pp)

Small and focused: 3–4 composite figures, each tied to a specific theorem, prose kept to one
paragraph per figure. Candidates (from Codex `07-numerics.tex`, re-voiced; protocols →
Appendix I):

1. **Three-stage trajectory** — burn-in linear in $\logp(1/(\beta\lambda_{\min}(C_0)))$
   across a scaling sweep + a hard $\kappa$-large run showing all three stages (ties to
   §2.2, §3.6); doubles as the paper's visual teaser if promoted to §1.
2. **Riemannian vs KL** — matched-step comparison + larger-step covariance robustness of the
   resolvent (ties to §2.3–2.4).
3. **STL floors** — paired-seed variance diagnostics and constant-step floors vs batch,
   floor vanishing for Gaussian targets (ties to §4.3–4.5).
4. **Non-log-concave boundary** — cascade, KL pole with the $4\varepsilon^{-2}$ reference,
   clipped-KL running-min under the $(h/N)$ envelope (ties to §2.7).

The affine-mode preconditioning experiment is kept only if space permits, as a single panel
attached to §5's modal-rate theorem (else cut).

### §7 Discussion (~1 p)

Two paragraphs of recap (what is sharp, what is only sufficient) + open problems: fixed-
dimensional local rate; sharpness of the $\sqrt{N_\theta}\log\kappa_\star$ branch; closed-form
discrete local radius; adaptive/implicit stepsizes vs the fixed-step lower bound; a stochastic
clipped stationarity theory (deliberately not claimed here); sharp rate theory for unbalanced
affine-invariant metrics; WFR dynamics as the natural second paper (one sentence, no
derivation).

---

## 3. Appendix architecture

| App | content | source |
|---|---|---|
| A | Gaussian differential identities; second-order score calculus ($\sqrt6$ chaos constant); $K_{ng}$ derivation | det/stoch note appendices |
| B | Logarithmic-spiral construction: Hessian certificate, shadowing, $\kappa_\star$ bounds, proof of the continuous lower bound | det note l.980–1715 (Codex `appendices/B` as condensed base) |
| C | Bump-train construction: flat-top bump, Gaussian averages, shadowing for both maps, proof of the discrete lower bound | det note l.2389–3010 |
| D | Convex-ridge local counterexample: tent profile, radial concentration, isotropic tuning, Codazzi conversion, Rayleigh direction, smooth realization | `local-log-kappa-counterexample.tex` (Codex `appendices/C` as condensed base) |
| E | Deterministic one-step proofs: scalar inequalities, retraction smoothness, KL/Bregman one-step (variational characterization, forward-KL descent) | det/stoch notes |
| F | Stochastic proofs: one-step lemmas, global theorems, matrix-exponential perturbation, local perturbation/killed-process/containment, decreasing-step | stoch note §12–15 + appendix |
| G | Non-log-concave details: double-well estimates, relative smoothness of the split divergence, clipping = log-det Bregman projection, mean-block identity for the $\lambda_+/\lambda_-$ factor | Codex `appendices/D` |
| H | Affine-invariant classification: irreducible $O(N_\theta)$-representation splitting, Schur argument | Codex `appendices/A` |

| I | Experimental protocols: targets, grids, seeds, quadrature, stopping rules (condensed from Codex `appendices/E`; keep the provenance-limits candor, condensed) | Codex `appendices/E` |

---

## 4. Source-of-truth rules and known defects to fix while writing

**Authority hierarchy (math):**
1. `improved-global-local.tex` — deterministic global/local rates, whitened/affine-invariant
   forms, both lower-bound constructions, exact Gaussian maps and gates.
2. `improved-global-local-stoch.tex` — everything stochastic (§10–17 + appendices), rescue,
   `lem:KL-logtangent`; its deterministic first half is a restatement — when the two notes
   differ, use the sharper/cleaner statement (the det note's whitened `thm:cont-global`
   supersedes the stoch note's original-coordinate version).
3. `local-log-kappa-counterexample.tex` — §3.5 + Appendix D, in the two-sided
   integration-ready form.
4. Codex draft — quarry only, three legitimate veins: (a) `06-nonconvex.tex` +
   `appendices/D` (double-well propositions; *repaired* clipped-KL theorem + $\lambda_+/\lambda_-$
   corollary — supersedes `natural-gradient.tex`'s clipped theorem); (b) `02-geometry.tex` +
   `appendices/A` (formal modal-rate + balanced-uniqueness theorems, Schur appendix);
   (c) condensed appendix write-ups of B/C/D as drafting bases. Everything else in the Codex
   draft is a re-voicing of the notes and is not used as text.
5. `natural-gradient.tex` — authoritative for *voice*, scheme derivations, the 1D example
   `rem:FR-init-cov`, the continuous stationarity corollary, the $(\eta,\omega,\tau)$ flow
   and its discretizations. Its global/local rate theorems and plain-estimator stochastic
   theorems are superseded and not carried.
6. `skeleton.tex` — fully superseded; not used.

**Defects to fix on import:**
- `natural-gradient.tex` l.155 red TODO tcolorbox (FR→W2 factor 1/2): mooted — replaced by
  the det note's corrected two basic inequalities.
- `natural-gradient.tex` l.829 STL $\hat K_n$ display typo.
- `natural-gradient.tex` clipped-KL theorem proof gap (clipping vs three-point lemma):
  superseded by the Bregman-projection proof (Codex Appendix D vein).
- `natural-gradient.tex` `rem:FR-init-cov` final sentence references the cut WFR section:
  rewrite to point at §2.6.
- ng `thm:affine-invariant-metrics` proof quotes the SPD classification with a duplicated
  $\Tr(C^{-1}XC^{-1}Y)$ term (l.2826): fix to $\Tr(C^{-1}X)\Tr(C^{-1}Y)$.
- Counterexample note's headline must be stated against the *new* spectral theorem
  (two-sided $\Theta$ form), i.e., un-comment and adopt the integration-ready corollary.

**Label/notation discipline:** one namespaced label scheme fixed before drafting
(e.g. `thm:glob-cont`, `thm:glob-riem`, `thm:loc-flow`, `thm:three-stage`,
`thm:stoch-*`, `thm:aff-*`); all imported labels renamed on entry; `\cref` throughout.

---

## 5. Resolved items (2026-07-22)

- **Q1 Numerics:** small focused section (§6) + Appendix I; Discussion becomes §7.
- **Q2 Stochastic clipped scheme:** dropped; open direction in §7.
- **Q3 Notation:** `natural-gradient.tex` conventions.
- **Q4 Drafting mechanics:** fresh rewrite of `sections/*.tex` in the author's voice; Codex
  files kept aside as quarry only (its three original veins: nonconvex §2.7+App G,
  geometry theorems §5+App H, condensed appendix bases B/C/D).

## 6. Drafting order (once the plan is approved)

1. `preamble.tex` refresh: notation macros aligned to `natural-gradient.tex`; label
   namespace fixed.
2. §2 (global) — largest section, anchors notation and the two schemes.
3. §3 (local) — including the counterexample statement; then §3.6 three-stage headline.
4. §4 (stochastic) — restating §3 content it uses, per [INSIGHT: S4↔S3].
5. §5 (geometry), §2.7 (nonconvex) — the two quarry-dependent units.
6. §1 (introduction) and abstract — written last, against the final numbered statements.
7. §6 numerics + §7 discussion; appendices in parallel with their sections.
8. Compile checks + label/citation audit per section, not at the end.
