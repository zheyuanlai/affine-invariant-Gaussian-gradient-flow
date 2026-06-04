# GPU backend on Colab Pro / A100 — operating notes

The PyTorch GPU backend reproduces the *corrected* CPU estimators (symmetrized
`H_sym`, Fisher–Rao-scaled `L_star`, diagonal benchmark) on a CUDA device so the
sample-size-scaling and operator/local-rate sweeps run at production scale. The
NumPy/SciPy CPU path is unchanged and remains the default; the torch path is
selected only when `operator.backend` resolves to `torch`.

> The torch backend is numerically identical to the NumPy path on the same bank
> (it copies `rho`, `M`, and the feature parameters from the CPU potential), so
> CPU and GPU results are directly comparable.

## Setup (Colab)

```bash
pip install -r requirements.txt
# Install a CUDA build of torch matching the Colab runtime (cu121 shown):
pip install torch --index-url https://download.pytorch.org/whl/cu121
```

### Device check

```bash
python -c "import torch; print(torch.cuda.is_available()); print(torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'no cuda')"
```

Or the bundled helper (prints the full device info, then runs the GPU smoke):

```bash
python scripts/natural_gradient_local_rate/colab_gpu_smoke.py --device cuda
```

## Recommended commands

**Smoke (fast, end-to-end sanity):**

```bash
python scripts/natural_gradient_local_rate/run_sample_size_scaling.py \
    --config configs/natural_gradient_local_rate/gpu_smoke.yaml \
    --backend torch --device cuda --overwrite
```

**Production sample-size scaling (run explicitly; this is the expensive sweep):**

```bash
python scripts/natural_gradient_local_rate/run_sample_size_scaling.py \
    --config configs/natural_gradient_local_rate/gpu_sample_size_scaling.yaml \
    --backend torch --device cuda --overwrite
```

The operator-grid and linearized-rate-grid runners take the same
`--backend/--device/--dtype/--chunk-size/--basis-block-size/--explicit-dense-max-N-theta`
flags (CLI overrides win over the YAML `operator:` block).

## How the GPU path works

* For `N_theta <= explicit_dense_max_N_theta` (default 64) the backend builds the
  dense `H_sym` (`p x p`, `p = N(N+1)/2`) and `L_star` (`D x D`, `D = N + p`)
  matrices as **chunked matmuls** over the sample bank — no Python basis loop —
  then calls `torch.linalg.eigh` entirely on device. `Lambda_hat_full_sym` is the
  top eigenvalue of `H_sym`; `gamma_loc` is the smallest eigenvalue of `L_star`.
* The diagonal benchmark `A = G - 11ᵀ` and the separable Gauss–Hermite exact
  benchmark are computed alongside (the exact benchmark reuses the cheap CPU
  quadrature).
* SciPy is never used on the GPU path, and there are no CPU/GPU transfers inside
  the eigensolve.

## Memory notes (float64 dense `L_star`, `D = N + N(N+1)/2`)

| `N_theta` | `p`   | `D`    | dense `L` memory | feasibility |
|----------:|------:|-------:|-----------------:|-------------|
| 16        | 136   | 152    | ~0.2 MB          | trivial     |
| 32        | 528   | 560    | ~2.5 MB          | trivial     |
| 64        | 2080  | 2144   | ~37 MB           | easy on A100 |
| 128       | 8256  | 8384   | ~560 MB          | feasible on A100; construction slower |

`explicit_dense_max_N_theta` defaults to 64. Raising it to 128 is feasible on an
A100 (the `eigh` workspace adds a few GB); beyond that, prefer the matrix-free
subspace iteration (present but not the supported production path in this pass).
Per-point peak GPU memory is recorded in the `gpu_peak_memory_mb` CSV column.

`device="cuda"` raises a clear error if CUDA is unavailable — it never silently
falls back to CPU. Use `device="cpu"` (with `backend="torch"`) to exercise the
torch path on a CPU-only machine, or `backend="auto"` to use torch only when a
CUDA device is present.

## Interpreting results (unchanged from the CPU analysis)

Do **not** read `Lambda_hat_full_sym` as a property of the flow until the
sample-size scaling has plateaued and the separable diagonal/exact benchmarks
agree. For a separable control the exact benchmark is `M_mc`-independent; if
`Lambda_hat_full_sym` keeps shrinking toward it as `M_mc` grows, the
high-dimensional inflation is finite-sample spectral noise, not real dimension
dependence. The GPU backend changes only *how fast* the estimates are computed,
not *what* they mean.
