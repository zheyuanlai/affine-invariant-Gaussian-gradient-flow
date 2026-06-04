import os
import subprocess
import sys

import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
SCRIPT = os.path.join(_ROOT, "scripts", "natural_gradient_local_rate",
                      "plot_estimator_diagnostics.py")


def test_diagnostic_plotting_smoke_no_ratio_by_default(tmp_path):
    base = tmp_path / "nglr"
    scaling = base / "sample_size_scaling"
    scaling.mkdir(parents=True)
    rows = []
    for M in [256, 1024]:
        for fam, lam, gamma, exact in [
            ("gaussian", 0.20 * (256 / M) ** 0.5, 0.8 + 0.15 * (M / 1024), 0.0),
            ("separable", 0.30 + 0.20 * (256 / M) ** 0.5, 0.75 + 0.10 * (M / 1024), 0.30),
            ("random_feature", 0.40 + 0.22 * (256 / M) ** 0.5, 0.70 + 0.12 * (M / 1024), float("nan")),
        ]:
            rows.append({
                "potential_family": fam,
                "family": fam,
                "N_theta": 4,
                "kappa_target": 5.0,
                "seed": 0,
                "M_mc": M,
                "Lambda_hat_full_sym": lam,
                "Lambda_hat_raw_forward": lam,
                "Lambda_hat_diag": lam / 2,
                "Lambda_hat_separable_exact": exact,
                "gamma_loc": gamma,
                "inverse_gamma_loc": 1.0 / gamma,
                "self_adjoint_error_H_sym": 1e-14,
                "self_adjoint_error_L_star": 1e-14,
                "full_over_diag": 2.0,
            })
    pd.DataFrame(rows).to_csv(scaling / "results_long.csv", index=False)

    outdir = tmp_path / "figures"
    res = subprocess.run(
        [sys.executable, SCRIPT, "--input", str(base), "--outdir", str(outdir)],
        cwd=_ROOT, capture_output=True, text=True, timeout=300,
    )
    assert res.returncode == 0, f"plot script failed:\nSTDOUT\n{res.stdout}\nSTDERR\n{res.stderr}"
    assert (outdir / "diag1_lambda_full_sym_vs_Mmc.png").exists()
    assert (outdir / "diag2_lambda_with_baseline_overlay.pdf").exists()
    assert (outdir / "diag10_gap_vs_sqrt_p_over_M.png").exists()
    assert not (outdir / "deprecated_full_over_diag_vs_Ntheta.png").exists()
