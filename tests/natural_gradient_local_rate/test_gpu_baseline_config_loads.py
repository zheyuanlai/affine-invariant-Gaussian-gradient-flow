from pathlib import Path

from src.common.io_utils import load_yaml


ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = ROOT / "configs" / "natural_gradient_local_rate"


def test_kept_nglr_configs_load_and_have_required_fields():
    """Every config under configs/natural_gradient_local_rate/ loads and is sane."""
    for path in sorted(CONFIG_DIR.glob("*.yaml")):
        cfg = load_yaml(path)
        assert cfg["experiment"]["name"], path.name
        assert cfg["grid"]["potential_family"], path.name
        assert "monte_carlo" in cfg, path.name
        assert cfg["outputs"]["base_dir"] == "outputs/natural_gradient_local_rate", path.name


def test_gpu_smoke_config_is_torch_dense():
    cfg = load_yaml(CONFIG_DIR / "gpu_smoke.yaml")
    op = cfg["operator"]
    assert op["backend"] == "torch"
    assert op["eigensolver"] == "torch_dense_eigh"
    assert op["explicit_dense_max_N_theta"] >= max(cfg["grid"]["N_theta"])
    assert isinstance(cfg["monte_carlo"]["M_mc"], list)


def test_lowdim_operator_full_config_is_the_final_production_scan():
    """The single production grid interpreted as evidence in the report."""
    cfg = load_yaml(CONFIG_DIR / "gpu_lowdim_operator_full.yaml")
    assert cfg["grid"]["N_theta"] == list(range(1, 17))
    assert cfg["grid"]["kappa_target"] == [2, 5, 10, 20, 50, 100]
    assert cfg["grid"]["seeds"] == [0, 1, 2]
    assert cfg["grid"]["potential_family"] == [
        "gaussian", "separable", "additive_index", "random_feature", "radial_tail",
    ]
    assert cfg["monte_carlo"]["M_mc"] == 4_194_304
    op = cfg["operator"]
    assert op["backend"] == "torch"
    assert op["eigensolver"] == "torch_dense_eigh"
    assert op["compute_gamma_loc"] is True
    assert op["explicit_dense_max_N_theta"] == 16
    assert op["chunk_size"] == 131_072
    assert cfg["monte_carlo"]["chunk_size"] == 131_072
    assert cfg["linearized_rate"]["save_eigenvectors"] is True
