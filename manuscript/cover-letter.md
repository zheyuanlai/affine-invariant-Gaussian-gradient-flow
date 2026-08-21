# Cover letter draft — Foundations of Computational Mathematics

*(Fill in author names, date, and the editor's name before submission. Verify
the live FoCM/EMS submission route and its current requirements immediately
before submitting; the journal is transitioning from Springer to EMS Press.)*

Dear Editors,

We submit "Gaussian Natural Gradient Flows for Variational Inference" for
consideration in Foundations of Computational Mathematics.

Gaussian variational inference by natural-gradient (Fisher–Rao) flow is a
widely used computational method whose convergence behavior has so far been
understood only through smoothness surrogates or particular covariance
parameterizations. Our paper gives what we believe is the first sharp,
intrinsic convergence theory for this method. The organizing result is that
convergence proceeds in three mathematically distinct stages — a covariance
burn-in logarithmic in the initialization, a global localization phase governed
by the condition number after whitening by the *optimizing* covariance, and a
local phase at the spectral gap of a self-adjoint score operator — and that
each stage's rate is sharp, certified by explicit constructions (a logarithmic
spiral, a bump train, and a convex-ridge family that rules out every
dimension-free logarithmic local rate). The same three-stage structure is
proved for two deterministic discretizations and for stochastic algorithms
built from a joint Price/Hessian oracle, with pathwise covariance bands and
explicit variance floors. The paper closes with a classification of all
affine-invariant metrics on the Gaussian manifold, identifying Fisher–Rao as
the unique member, up to scale, balancing the mean, trace, and traceless
covariance modes.

We believe the paper fits FoCM's remit — the mathematics underlying
computation — because it connects information geometry, matrix dynamics,
sharp complexity lower bounds, and stochastic approximation into a single
quantitative account of a practical algorithm.

The manuscript is not under consideration elsewhere.

**AI-use disclosure** (per EMS policy; adjust wording as appropriate): AI
assistants (Anthropic Claude, OpenAI Codex/ChatGPT) were used to assist with
manuscript drafting, consistency auditing, and software development. The
authors have reviewed and verified all mathematical statements, proofs,
computations, citations, and numerical results, and take full responsibility
for the content.

Sincerely,
[Corresponding author, affiliation, email]
