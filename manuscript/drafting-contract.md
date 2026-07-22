# Drafting contract — binding on every section draft

Read this before writing any LaTeX. Companion to `chapter-plan.md` (structure) and
`natural-gradient.tex` (voice).

---

## 1. Voice — imitate `natural-gradient.tex`, the author's own writing

The paper must read as though written by the author of `natural-gradient.tex`. Concretely:

**Do:**
- First-person plural, walking narrative. One lead-in sentence before each result:
  > "We proceed to show the convergence theory of \eqref{upd:Riem} under Assumption \ref{assump:logconcave-smooth}."
  > "We first show the spectral bounds of the covariance along \eqref{ODEs:NG}."
  > "We then prove the spectral bounds of covariance along \eqref{upd:Riem}."
- One pragmatic sentence of motivation, then the mathematics:
  > "The discrete algorithm with Riemannian distance \eqref{upd:Riem} is geometrically natural but computationally expensive, which requires to compute matrix exponential at every step, for the ease in computation, we develop the following discretization scheme with KL divergence."
- Hypotheses inline in the statement: "Under Assumption \ref{assump:logconcave-smooth}, if the
  initial covariance satisfies $\lambda_{0,\min} I \preceq C_0 \preceq \lambda_{0,\max} I$, then …"
- Proofs as direct computation with short connectives: "This gives", "it follows that",
  "by Gronwall's inequality", "which completes the proof".
- Derivations shown in the flow of the text (forward step, backward step, then the boxed
  update), as in `natural-gradient.tex` §2.1–2.2.

**Do not:**
- No bold run-in paragraph headers (`\paragraph{A covariance bootstrap...}`) in the body.
  That is the *notes'* voice, not the paper's.
- No meta-commentary about the document: no "this note records", "we now turn to", "in this
  section we will discuss", "the reader should note", "it is important to note".
- No proof-strategy previews ("The proof proceeds in three steps: first we…"). Just prove it.
- No defensive hedging blocks. Scope limitations become **one short `remark`** after the
  theorem, in plain declarative sentences.
- No em-dash-heavy prose; the author uses commas and short clauses.
- No bullet lists in the body except where genuinely enumerating cases of a definition.

**Register calibration.** Read `natural-gradient.tex` lines 87–120, 306–400, 504–545 before
drafting. Those three passages fix the voice.

---

## 2. Notation — `natural-gradient.tex` conventions, no exceptions

| object | symbol | macro |
|---|---|---|
| dimension | $N_\theta$ | `N_\theta` |
| variable | $\theta$ | |
| target | $\rho_\post(\theta) \propto \exp(-V(\theta))$ | `\rho_\post` |
| Gaussian family | $\PG$ | `\PG` |
| parameter | $a=(m,C) \in \Aspace = \R^{N_\theta}\times\SPD(N_\theta)$ | `\Aspace` |
| variational density | $\rho_a = \N(m,C)$ | |
| energy | $\calE(a) = \KL(\rho_a\|\rho_\post)$ | `\calE` |
| optimizer | $a_\star=(m_\star,C_\star)$ | |
| energy gap | $\DeltaE(a) = \calE(a)-\calE(a_\star)$, $\DeltaE_0 = \DeltaE(a_0)$ | `\DeltaE` |
| stepsize | $\Delta t$ | |
| iterate index | $n$; horizon $N$ | |
| curvature bounds | $\alpha \Id \preceq \nabla^2 V \preceq \beta \Id$, $\kappa=\beta/\alpha$ | |
| covariance bounds | $\lambda_{\min},\lambda_{\max}$; initial $\lambda_{0,\min},\lambda_{0,\max}$ | |
| whitened curvature | $\alphastar,\betastar,\kappastar$ | `\alphastar` etc. |
| whitened initial covariance | $\lamzstar = \lambda_{\min}(C_\star^{-1/2}C_0C_\star^{-1/2})$ | `\lamzstar` |
| local spectral scale | $\Gamma$, gap $\gstar$, top $\Lstar$ | `\gstar`, `\Lstar` |
| Fisher–Rao gradient | $\grad \calE(a)$, norm $\norm{\grad\calE(a)}_a$ | |

**Translation rules when importing from the notes** (they use $d$, $\pi$, $h$):
`d` → `N_\theta`, `\pi` → `\rho_\post`, `h` → `\Delta t`, `x` → `\theta` (integration variable),
`\Delta(a)` stays `\DeltaE(a)`. Keep $\Gamma$, $\alphastar$, $\betastar$, $\kappastar$,
$\lamzstar$, $\LHs$, $\barLH$, $\Kng$, $\Kdisc$ exactly as in the notes.

Macros are defined in `preamble.tex` — use them; do not redefine anything locally.

---

## 3. Label namespace — fixed, use exactly these

Rename every imported label on entry. Pattern: `type:section-slug`.

**Assumptions.** `assump:logconcave-smooth` (α,β), `assump:smooth` (bounded Hessian only,
non-log-concave), `assump:hess-lip` (Hessian-Lipschitz).

**Equations (schemes and flows).** `ODEs:NG`, `upd:Riem`, `upd:KL`, `upd:clip-KL`,
`upd:Riem-stoch`, `upd:KL-stoch`, `ODEs:affine`, `upd:affine-Riem`, `upd:affine-KL`.

**§2 global.** `lem:cov-cont`, `thm:glob-cont`, `rem:FR-init-cov`, `thm:sharp-cont`,
`prop:upd-map-Riem`, `lem:cov-Riem`, `lem:descent-Riem`, `thm:glob-Riem`, `lem:cov-KL`,
`prop:smooth-Bregman`, `prop:onestep-KL`, `thm:glob-KL`, `thm:sharp-disc`, `cor:sharp-disc`,
`lem:cov-BW`, `thm:burnin-BW`, `def:rescue`, `lem:rescue-band`, `thm:rescue-Gauss`,
`thm:burnin-rescue`, `cor:stationary-cont`, `prop:cascade-Riem`, `prop:pole-KL`,
`thm:clip-KL`, `cor:clip-KL-full`.

**§3 local.** `def:whitened`, `def:score-ops`, `lem:diagonal-modes`, `prop:spectral-sandwich`,
`prop:linearized-maps`, `lem:coercivity`, `def:region`, `lem:gauss-core`, `lem:Kng`,
`thm:loc-flow`, `cor:loc-universal`, `cor:loc-Gauss`, `prop:Gauss-KL-map`,
`prop:Gauss-Riem-map`, `thm:loc-disc`, `lem:loc-contraction`, `prop:entry-residual`,
`prop:entry-Gauss`, `thm:loc-sharp`, `thm:three-stage`, `cor:three-stage-burninfree`.

**§4 stochastic.** `def:oracle`, `alg:unified`, `lem:stoch-band`, `def:Psi`, `lem:Psi-selfbound`,
`lem:Psi-hesslip`, `lem:onestep-stoch-Riem`, `lem:onestep-stoch-KL`, `thm:stoch-glob-Riem`,
`thm:stoch-glob-KL`, `lem:KL-retraction`, `lem:loc-variance`, `lem:perturb-Riem`,
`lem:perturb-KL`, `thm:stoch-loc-Riem`, `thm:stoch-loc-KL`, `prop:exit-prob`,
`rem:stoch-loc-complexity`, `cor:stoch-loc-Gauss`, `cor:stoch-three-stage`,
`thm:decreasing`, `prop:Gauss-nofloor`, `tab:complexity`.

**§5 geometry.** `thm:aff-class`, `rem:aff-dim-one`, `thm:aff-flow`, `thm:aff-modes`,
`thm:aff-balanced`, `rem:aff-scope`.

**Sections.** `sec:intro`, `sec:global`, `sec:local`, `sec:stochastic`, `sec:geometry`,
`sec:numerics`, `sec:discussion`. **Appendices.** `app:identities`, `app:spiral`,
`app:bump`, `app:ridge`, `app:det-proofs`, `app:stoch-proofs`, `app:nonconvex`,
`app:classification`, `app:experiments`.

Cross-reference with `\Cref{...}` at sentence start, `\cref{...}` inside a sentence.
Equations with `\eqref{...}`.

---

## 4. Mathematical source of truth

Authority order — when two sources disagree, the higher one wins:

1. **`improved-global-local.tex`** — deterministic global/local rates, whitened forms, both
   lower-bound constructions, exact Gaussian maps and entry gates.
2. **`improved-global-local-stoch.tex`** — everything stochastic (its §10–17 and appendices),
   the quadratic rescue, and `lem:KL-logtangent`. Its *first half* is a restatement of source 1;
   when both state a result, take the sharper/cleaner one (e.g. the whitened `thm:cont-global`
   of source 1 supersedes the original-coordinate version here).
3. **`local-log-kappa-counterexample.tex`** — §3.5 and Appendix D, stated in the two-sided
   integration-ready form (its commented-out corollary at lines 1661–1670).
4. **Codex draft** (`sections/`, `appendices/`) — quarry, three veins only:
   `06-nonconvex.tex` + `appendices/D` (double-well propositions; repaired clipped-KL theorem
   and its $\lambda_+/\lambda_-$ corollary); `02-geometry.tex` + `appendices/A` (formal modal-rate
   and balanced-uniqueness theorems, Schur classification); condensed appendix write-ups as
   drafting bases. Its *prose is not reused* — re-voice everything.
5. **`natural-gradient.tex`** — authoritative for voice, the two scheme derivations, the 1D
   example `rem:FR-init-cov`, the continuous stationarity corollary, and the $(\eta,\omega,\tau)$
   flow with its discretizations. Its global/local rate theorems and plain single-sample
   stochastic theorems are **superseded — do not carry them**.
6. `skeleton.tex` — superseded, unused.

**IRON RULE — no invented mathematics.** Every theorem, constant, and rate must be traceable
to one of the sources above. If a needed statement does not exist in any source, do not
invent it: insert
`% TODO(gap): <precise description of the missing statement and why it is needed>`
and continue. Do not silently weaken, strengthen, or "clean up" a constant.

**Defects to fix on import** (already diagnosed):
- The FR→$W_2$ gradient-norm comparison in `natural-gradient.tex` (red TODO tcolorbox, l.155)
  is superseded by the corrected two basic inequalities in source 1. Drop the tcolorbox.
- `natural-gradient.tex` l.829: the STL $\hat K_n$ display has a typo; take the estimator
  definition from the stochastic note instead.
- `natural-gradient.tex`'s clipped-KL theorem has a proof gap (clipping vs the three-point
  lemma); use the Codex/Bregman-projection version.
- `natural-gradient.tex` l.2826: the quoted SPD classification duplicates
  $\Tr(C^{-1}XC^{-1}Y)$; the second term is $\tau\Tr(C^{-1}X)\Tr(C^{-1}Y)$.
- `rem:FR-init-cov`'s final sentence points at the cut WFR section; rewrite it to point at the
  burn-in removal subsection (§2.6).

---

## 5. Statement economy

Main text carries **headline statements plus what is needed to state them**. Everything else
goes to an appendix, invoked by one sentence.

- A section should have on the order of 5–10 numbered main-text statements, not 30.
- Long or technical proofs (>~1 page): state in main text, prove in the appendix with
  "The proof is given in \Cref{app:...}."
- Short proofs (<~half a page) stay in place — the author's manuscript proves in place.
- Supporting lemmas that only serve one proof belong in that proof's appendix.
- Complexity statements go **inside** the theorem: clean rate first, then "In particular, with
  $\Delta t \asymp \ldots$, we obtain $\DeltaE(a_N) \leq \varepsilon$ after
  $N = O(\ldots)$ iterations."

---

## 6. Cross-section referencing (author's instruction, load-bearing)

Section 4 must **restate the content** of the Section 3 results it uses, not merely name them.
Write, for example:

> "The one-step contraction of \Cref{lem:loc-contraction} holds pointwise on $\Ureg_\delta$:
> every $a \in \Ureg_\delta$ satisfies $\norm{T(a)-a_\star} \le (1-\tfrac12\Delta t\,\gamma_\bullet)\norm{a-a_\star}$,
> and requires only $\Delta t \le 2/(4+\Gamma)$, with no energy-descent restriction on the
> stepsize. We perturb this estimate by the Price/Hessian noise."

not

> "By \Cref{lem:loc-contraction} and \Cref{prop:spectral-sandwich}, we obtain …"

The same rule applies to any backward reference carrying a rate, constant, radius, or
stepsize window.

---

## 7. Citations

`refs.bib` is the only bibliography. Cite with `\citet`/`\citep` (natbib, numbered).
Key entries: `lambert2022variational`, `diao2023forward`, `carrillo2026fisher`,
`chen2023sampling`, `amari1998natural`, `amari2016information`, `rao1945information`,
`blei2017variational`, `JMLR:v14:hoffman13a`, `roeder2017sticking`, `thanwerdas2019affine`,
`thanwerdas2023n`, `sun2025natural`, `parikh2014proximal`, `brooks2011handbook`.
**Never invent a citation key.** If a needed reference is absent from `refs.bib`, write
`% TODO(cite): <what is needed>`.

---

## 8. Output mechanics

- One file per section, `\section{...}\label{sec:...}` at the top, no preamble, no
  `\begin{document}`.
- Appendix files: `\section{...}\label{app:...}` (the `\appendix` switch lives in `main.tex`).
- The file must compile as part of `main.tex` — balanced environments, all macros from
  `preamble.tex`.
- End every drafted file with a comment block listing the sources used, in the form
  `% SOURCES: improved-global-local.tex l.885 (thm:cont-global) -> thm:glob-cont`, so the
  provenance of each statement is auditable.
