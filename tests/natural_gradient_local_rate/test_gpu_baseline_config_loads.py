from pathlib import Path

from src.common.io_utils import load_yaml


ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = ROOT / "configs" / "natural_gradient_local_rate"


def test_gpu_baseline_configs_load_and_include_required_fields():
    for name in [
        "gpu_baseline_scaling.yaml",
        "gpu_N64_deep_scaling.yaml",
        "gpu_lowdim_highM_scaling.yaml",
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


def test_lowdim_highM_stage1_config_is_high_sample_calibration():
    cfg = load_yaml(CONFIG_DIR / "gpu_lowdim_highM_scaling.yaml")
    assert cfg["grid"]["N_theta"] == [8, 16, 24, 32]
    assert cfg["grid"]["kappa_target"] == [5]
    assert cfg["grid"]["seeds"] == [0, 1, 2]
    assert cfg["grid"]["potential_family"] == [
        "gaussian", "separable", "random_feature", "radial_tail",
    ]
    assert max(cfg["monte_carlo"]["M_mc"]) >= 4_194_304
    assert min(cfg["monte_carlo"]["M_mc"]) >= 262_144
    assert cfg["operator"]["explicit_dense_max_N_theta"] == 32
    assert cfg["operator"]["chunk_size"] == 8192
