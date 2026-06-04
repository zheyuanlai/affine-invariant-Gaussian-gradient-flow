from pathlib import Path

from src.common.io_utils import load_yaml


ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = ROOT / "configs" / "natural_gradient_local_rate"


def test_gpu_baseline_configs_load_and_include_required_fields():
    for name in [
        "gpu_baseline_scaling.yaml",
        "gpu_N64_deep_scaling.yaml",
        "gpu_smoke_gaussian.yaml",
    ]:
        cfg = load_yaml(CONFIG_DIR / name)
        assert cfg["experiment"]["name"]
        assert "gaussian" in cfg["grid"]["potential_family"]
        assert "separable" in cfg["grid"]["potential_family"]
        assert "random_feature" in cfg["grid"]["potential_family"]
        assert isinstance(cfg["monte_carlo"]["M_mc"], list)
        op = cfg["operator"]
        assert op["backend"] == "torch"
        assert op["eigensolver"] == "torch_dense_eigh"
        assert op["compute_gamma_loc"] is True
        assert op["explicit_dense_max_N_theta"] >= max(cfg["grid"]["N_theta"])
        assert cfg["outputs"]["base_dir"] == "outputs/natural_gradient_local_rate"
