# Handoff: Logarithmic Dimension-Free Local Rate

## Current Objective

Work toward proving or refuting the logarithmic dimension-free local-rate conjecture for the equilibrium-whitened Gaussian natural-gradient flow:

```tex
E[\nabla V(Z)] = 0,
E[\nabla^2 V(Z)] = I,
\alpha I \preceq \nabla^2 V(Z) \preceq \beta I,
\kappa=\beta/\alpha,
```

with desired rate

```tex
-\lambda_{\star,\max} \gtrsim 1/(1+\log \kappa)
```

independent of `N_theta`.

## Mathematical State

The old route `Lambda <= D_kappa - 1` should not be pursued. The radial threshold construction in `local-convergence/fable-report.tex` gives a counterexample to that fourth-order operator-norm bound, but that example is even/radial, has `T=0`, and does not slow the actual local rate.

The useful new route is through the first-Hermite Hessian coupling:

```tex
\tau_H := \|T\|_{\mathrm{op}}
       = \sup_{\|w\|=1} \| E[\nabla^2 V(Z)(w^T Z)] \|_F.
```

`local-convergence/loc-conv-sec.tex` now contains a deterministic reduction:

```tex
\tau_H^2 <= M  ==>  -\lambda_{\star,\max} >= 1/(M+3).
```

It also records the baseline dimension-free bound:

```tex
\tau_H^2 <= \beta,
\qquad
-\lambda_{\star,\max} >= 1/(\beta+3) >= 1/(\kappa+3).
```

This settles the weak dimension-free claim. It does not settle the desired logarithmic dependence.

So the remaining proof target is:

```tex
\tau_H^2 <= C(1+\log\kappa)
```

dimension-free in `N_theta`.

The hard analytic subproblem has been isolated as the transverse term, for fixed unit `w` and `Z = t w + y`:

```tex
sup_{||B||_F=1, Bw=0} | E[(B y)^T \nabla^2 V(Z) w] |.
```

The manuscript now proves the easy rank-one part:

```tex
|E[(w^T Z) q^T \nabla^2 V(Z) q]|
  <= C_0 sqrt(1+log kappa),
```

and by polarization

```tex
for unit w:
|E[(w^T Z) p^T \nabla^2 V(Z) q]|
  <= C_0 ||p|| ||q|| sqrt(1+log kappa).
```

Therefore the longitudinal and mixed pieces in the directional decomposition are already logarithmic. The unresolved issue is genuinely high-rank: the transverse trace/contraction against `B`.

A counterexample must make this true Hessian-compatible coupling large and must show `gamma_loc * (1+log kappa) -> 0`. A large positive covariance `Lambda` is not enough.

Important relaxed-field warning: if Hessian compatibility is dropped, `tau_H` can grow with dimension. The field

```tex
W(z)=I+\epsilon \tanh(z_1) I
```

has `E[W]=I` and bounded eigenvalues for small `epsilon`, but `||E[W z_1]||_F ~ epsilon sqrt(N_theta)`. It is not a Hessian because compatibility would require `partial_1 W_22 = partial_2 W_21`, while the left side is nonzero and the right side is zero. So the desired proof must use integrability of `W = Hess V`, not just ellipticity/whitening.

## Code Changes Already Made

New diagnostics:

- `src/natural_gradient_local_rate/operators.py`
  - Added dense `T_matrix(...)`.
  - Added batched Frobenius-isometric symmetric vectorization.
- `src/natural_gradient_local_rate/linearized_rate.py`
  - Added `estimate_tau_H(...)`.
  - Added top-`T` singular-mode decomposition into longitudinal/mixed/transverse blocks.
- `src/natural_gradient_local_rate/estimator_suite.py`
  - Rows now report:
    - `tau_H`
    - `tau_H_sq`
    - `coupling_bound_rate = 1/(tau_H_sq+3)`
    - `gamma_over_coupling_bound`
    - `tau_top_longitudinal`, `tau_top_mixed`, `tau_top_transverse`
    - corresponding fractions and top-mode block norm squares
- `src/natural_gradient_local_rate/torch_backend.py`
  - Torch backend now computes `tau_H` from the already-built `T_mat`.
  - Torch backend now reports the same top-`T` longitudinal/mixed/transverse decomposition.
  - Added torch support for `product_feature`.

New adversarial Hessian-compatible family:

- `src/natural_gradient_local_rate/potentials/product_feature.py`
  - Scalar feature:
    ```tex
    Phi(theta) = (1/sqrt(r)) sum_l a_l tanh(w_l^T theta+c_l) tanh(v_l^T theta+d_l)
    ```
  - Pairwise transverse directions `w_l`, `v_l`.
  - Intended to stress the transverse coupling while remaining a true scalar potential.

New config:

- `configs/natural_gradient_local_rate/adversarial_transverse.yaml`
  - CPU/default config for product-feature/random-feature/additive-index/radial-tail stress testing.
  - Override `--backend torch --device cuda` on GPU.

Runner sharding:

- `scripts/natural_gradient_local_rate/_common.py`
  - Added `--num-shards` and `--shard-index`.
- The operator, linearized-rate, and joint runners now write shard-suffixed CSVs, e.g.
  - `results_long_shard03-of-08.csv`
  - `summary_shard03-of-08.csv`

## Important Workspace Note

Git currently reports `local-convergence/` as untracked even though it contains the manuscript files. The edited manuscript file is:

```text
local-convergence/loc-conv-sec.tex
```

Make sure the GPU-side worktree includes this directory if manuscript edits matter.

## Processes

The long local CPU adversarial scan was stopped on request.

It had been writing partial files under:

```text
/tmp/nglr_adversarial_transverse
```

Do not treat those partial outputs as evidence.

## Verification Already Run

Before this handoff:

```bash
pytest tests/natural_gradient_local_rate -q
# 211 passed, 5 skipped
```

After adding torch support for `product_feature`:

```bash
pytest \
  tests/natural_gradient_local_rate/test_torch_operator_cpu_consistency.py \
  tests/natural_gradient_local_rate/test_torch_fast_paths.py \
  tests/natural_gradient_local_rate/test_tau_coupling_bound.py \
  -q
# 62 passed, 5 skipped
```

Shard smoke:

```bash
python scripts/natural_gradient_local_rate/run_operator_linearized_grid.py \
  --config configs/natural_gradient_local_rate/smoke.yaml \
  --outdir /tmp/nglr_shard_smoke \
  --overwrite \
  --num-shards 2 \
  --shard-index 1
```

This succeeded and wrote shard-suffixed CSVs.

After adding the rank-one Hermite lemma and the top-`T` decomposition diagnostics:

```bash
pytest \
  tests/natural_gradient_local_rate/test_tau_coupling_bound.py \
  tests/natural_gradient_local_rate/test_torch_operator_cpu_consistency.py \
  -q
# 34 passed, 5 skipped

pytest \
  tests/natural_gradient_local_rate/test_sample_size_scaling_smoke.py \
  tests/natural_gradient_local_rate/test_torch_runner_smokes.py \
  -q
# 5 passed

python -m compileall -q \
  src/natural_gradient_local_rate/linearized_rate.py \
  src/natural_gradient_local_rate/estimator_suite.py \
  src/natural_gradient_local_rate/torch_backend.py \
  scripts/natural_gradient_local_rate/_common.py
```

Runner smoke:

```bash
python scripts/natural_gradient_local_rate/run_operator_linearized_grid.py \
  --config configs/natural_gradient_local_rate/smoke.yaml \
  --outdir /tmp/nglr_tau_decomp_smoke \
  --overwrite
```

This completed all 24 rows. The new decomposition columns satisfy, in the smoke CSV, contribution additivity to `2.1e-16` and block-norm additivity to `1.8e-15`.

Small local CPU product-feature diagnostic:

```text
N in {4,8,12}, kappa in {10,100}, seeds {0,1}
```

No counterexample trend appeared. `tau_H_sq/(1+log kappa)` stayed around `0.0015-0.003`, while `Lambda_hat_full_sym` increased with dimension. This supports the distinction that large covariance curvature is not the slow-rate mechanism; it is not theorem-level evidence.

## Recommended H200 Run

On an 8-H200 machine, first verify the torch path:

```bash
pytest \
  tests/natural_gradient_local_rate/test_torch_operator_cpu_consistency.py \
  tests/natural_gradient_local_rate/test_tau_coupling_bound.py \
  -q
```

Then launch one process per GPU. From the repo root:

```bash
OUT=outputs/natural_gradient_local_rate/adversarial_h200
CFG=configs/natural_gradient_local_rate/adversarial_transverse.yaml
mkdir -p "$OUT"

for i in 0 1 2 3 4 5 6 7; do
  CUDA_VISIBLE_DEVICES=$i python scripts/natural_gradient_local_rate/run_operator_linearized_grid.py \
    --config "$CFG" \
    --outdir "$OUT" \
    --overwrite \
    --backend torch \
    --device cuda \
    --dtype float64 \
    --chunk-size 262144 \
    --explicit-dense-max-N-theta 64 \
    --num-shards 8 \
    --shard-index $i \
    > "$OUT/shard_$i.log" 2>&1 &
done
wait
```

If memory is tight, lower `--chunk-size` to `131072`. If utilization is low and memory is comfortable, try `524288`.

Combine shards:

```bash
python - <<'PY'
from pathlib import Path
import pandas as pd

base = Path("outputs/natural_gradient_local_rate/adversarial_h200")
lr = base / "linearized_rate_grid"
op = base / "operator_grid"

for folder in [lr, op]:
    frames = [pd.read_csv(p) for p in sorted(folder.glob("results_long_shard*.csv"))]
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(folder / "results_long.csv", index=False)
    print(folder, df.shape)
PY
```

## First Analysis to Run on GPU Outputs

```bash
python - <<'PY'
import numpy as np
import pandas as pd

p = "outputs/natural_gradient_local_rate/adversarial_h200/linearized_rate_grid/results_long.csv"
df = pd.read_csv(p)
df = df[df.status == "ok"].copy()
df["logfac"] = 1.0 + np.log(df["kappa_target"])
df["scaled_rate"] = df["gamma_loc"] * df["logfac"]
df["tau_sq_over_log"] = df["tau_H_sq"] / df["logfac"]

cols = [
    "gamma_loc", "scaled_rate", "tau_H_sq", "tau_sq_over_log",
    "coupling_bound_rate", "gamma_over_coupling_bound",
    "tau_top_transverse_fraction", "tau_top_X_transverse_norm_sq",
    "Lambda_hat_full_sym",
]
summary = df.groupby(["potential_family", "N_theta", "kappa_target"])[cols].agg(["mean", "min", "max"])
print(summary.to_string(max_rows=200))

print("\nWorst scaled_rate rows:")
print(df.sort_values("scaled_rate").head(20)[[
    "potential_family", "N_theta", "kappa_target", "seed",
    "gamma_loc", "scaled_rate", "tau_H_sq", "tau_sq_over_log",
    "gamma_over_coupling_bound",
    "tau_top_longitudinal_fraction", "tau_top_mixed_fraction",
    "tau_top_transverse_fraction", "tau_top_X_transverse_norm_sq",
    "Lambda_hat_full_sym",
    "empirical_min_hess_eig", "empirical_max_hess_eig",
]].to_string(index=False))
PY
```

Interpretation:

- If `scaled_rate = gamma_loc * (1+log kappa)` stays bounded below and `tau_sq_over_log` stays bounded, this supports the `tau_H` proof route.
- If `tau_sq_over_log` grows with `N_theta` for `product_feature`, inspect whether `gamma_loc` actually decreases. Large `tau_H` only matters if the slow eigenmode follows it.
- If `tau_top_transverse_fraction` is small, the run is not stressing the remaining hard term; it is mostly in the already-controlled rank-one/mixed part.
- If `tau_top_transverse_fraction` is large and `scaled_rate` decays, inspect that row first as the most plausible counterexample signal.
- A real counterexample needs a trend toward `scaled_rate -> 0`, not just growth in `Lambda_hat_full_sym`.

## Good Next Codex Tasks

1. Finish/verify multi-GPU adversarial scan on H200.
2. Inspect worst slow modes from saved eigenvectors, especially product-feature rows.
3. If no counterexample appears, update `reports/natural_gradient_local_rate_report.tex` with a new `tau_H` table/figure.
4. If `tau_H` grows, derive the corresponding analytic ansatz from `product_feature` and check whether it survives as `M_mc`, `N_theta`, and `kappa` increase.
5. Keep the mathematical target honest: either prove `tau_H^2 <= C(1+log kappa)` or find a Hessian-compatible sequence with `gamma_loc(1+log kappa) -> 0`.
