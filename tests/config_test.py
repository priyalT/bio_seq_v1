import pytest
import yaml
from pathlib import Path
from bio_seq_v1.config import Config


def test_config_loads_defaults():
    """Config should initialize with default values."""
    cfg = Config()
    assert cfg.config is not None
    assert isinstance(cfg.config, dict)


def test_config_source_is_defaults():
    """Without a config file present, source should be 'defaults'."""
    cfg = Config()
    assert cfg._config_source in ("defaults",) or Path(cfg._config_source).exists()


def test_defaults_has_all_sections():
    """Default config should contain all expected top-level sections."""
    cfg = Config()
    expected_sections = [
        "parsing",
        "output",
        "export",
        "statistics",
        "translation",
        "orf",
        "motif",
        "performance",
        "cli",
        "files",
        "logging",
    ]
    for section in expected_sections:
        assert section in cfg.config, f"Missing section: {section}"


def test_defaults_parsing_values():
    cfg = Config()
    parsing = cfg.config["parsing"]
    assert parsing["strict_mode"] is False
    assert parsing["strict_file"] is False
    assert parsing["strict_sequence"] is False
    assert parsing["max_errors"] == 100
    assert parsing["case_sensitive"] is False


def test_defaults_export_values():
    cfg = Config()
    export = cfg.config["export"]
    assert export["default_format"] == "csv"
    assert export["include_headers"] is True
    assert export["pretty_json"] is True
    assert export["json_indent"] == 2


def test_defaults_translation_values():
    cfg = Config()
    translation = cfg.config["translation"]
    assert translation["start_codons"] == ["ATG"]
    assert translation["stop_codons"] == ["TAA", "TAG", "TGA"]
    assert translation["unknown_codon_char"] == "X"


def test_defaults_orf_values():
    cfg = Config()
    orf = cfg.config["orf"]
    assert orf["min_length"] == 0
    assert orf["search_both_strands"] is True
    assert orf["nested_orfs"] is True


def test_defaults_motif_values():
    cfg = Config()
    motif = cfg.config["motif"]
    assert motif["default_k"] == 6
    assert motif["max_mismatches"] == 0
    assert motif["iupac_support"] is True



def test_repr():
    cfg = Config()
    r = repr(cfg)
    assert "Config(loaded_from=" in r



def test_get_config_simple_key():
    cfg = Config()
    value = cfg.get_config("export.default_format")
    assert value == "csv"


def test_get_config_nested_key():
    cfg = Config()
    value = cfg.get_config("output.html_report.template")
    assert value == "default"


def test_get_config_nonexistent_key_returns_default():
    cfg = Config()
    value = cfg.get_config("nonexistent.key", default="fallback")
    assert value == "fallback"


def test_get_config_nonexistent_key_returns_none():
    cfg = Config()
    value = cfg.get_config("nonexistent.key")
    assert value is None


def test_get_config_partial_path():
    """Getting a section path should return the whole sub-dict."""
    cfg = Config()
    value = cfg.get_config("parsing")
    assert isinstance(value, dict)
    assert "strict_mode" in value


def test_get_config_deep_nonexistent():
    cfg = Config()
    value = cfg.get_config("parsing.nonexistent.deep.key", default="nope")
    assert value == "nope"



def test_set_config_existing_key():
    cfg = Config()
    cfg.set_config("export.default_format", "json")
    assert cfg.config["export"]["default_format"] == "json"


def test_set_config_nested_key():
    cfg = Config()
    cfg.set_config("motif.default_k", 10)
    assert cfg.config["motif"]["default_k"] == 10


def test_set_config_creates_new_key():
    cfg = Config()
    cfg.set_config("custom.new_setting", "hello")
    assert cfg.config["custom"]["new_setting"] == "hello"


def test_set_config_creates_nested_path():
    cfg = Config()
    cfg.set_config("a.b.c.d", 42)
    assert cfg.config["a"]["b"]["c"]["d"] == 42


def test_set_config_overwrites_value():
    cfg = Config()
    cfg.set_config("parsing.max_errors", 50)
    assert cfg.config["parsing"]["max_errors"] == 50
    cfg.set_config("parsing.max_errors", 200)
    assert cfg.config["parsing"]["max_errors"] == 200


def test_deep_merge_simple():
    cfg = Config()
    base = {"a": 1, "b": 2}
    override = {"b": 3, "c": 4}
    cfg._deep_merge(base, override)
    assert base == {"a": 1, "b": 3, "c": 4}


def test_deep_merge_nested():
    cfg = Config()
    base = {"section": {"key1": "val1", "key2": "val2"}}
    override = {"section": {"key2": "updated", "key3": "new"}}
    cfg._deep_merge(base, override)
    assert base["section"]["key1"] == "val1"
    assert base["section"]["key2"] == "updated"
    assert base["section"]["key3"] == "new"


def test_deep_merge_override_replaces_non_dict():
    cfg = Config()
    base = {"key": "string_value"}
    override = {"key": {"nested": True}}
    cfg._deep_merge(base, override)
    assert base["key"] == {"nested": True}


def test_load_config_from_file(tmp_path):
    custom_config = {
        "parsing": {"strict_mode": True},
        "motif": {"default_k": 12},
    }
    config_file = tmp_path / "test_config.yaml"
    with open(config_file, "w") as f:
        yaml.dump(custom_config, f)

    cfg = Config()
    cfg.load_config(str(config_file))

    assert cfg.config["parsing"]["strict_mode"] is True
    assert cfg.config["motif"]["default_k"] == 12
    assert cfg.config["export"]["default_format"] == "csv"


def test_load_config_nonexistent_file(capsys):
    cfg = Config()
    cfg.load_config("/fake/nonexistent/config.yaml")
    captured = capsys.readouterr()
    assert "not found" in captured.out.lower()


def test_load_config_invalid_yaml(tmp_path, capsys):
    bad_yaml = tmp_path / "bad.yaml"
    bad_yaml.write_text(":\n  - :\n    invalid: [yaml: {{broken")

    cfg = Config()
    cfg.load_config(str(bad_yaml))
    captured = capsys.readouterr()
    assert "error" in captured.out.lower() or "yaml" in captured.out.lower()


def test_load_config_empty_file(tmp_path):
    empty_file = tmp_path / "empty.yaml"
    empty_file.write_text("")

    cfg = Config()
    original_config = cfg.config.copy()
    cfg.load_config(str(empty_file))
    assert cfg.config["parsing"]["strict_mode"] == original_config["parsing"]["strict_mode"]


def test_load_config_merges_not_replaces(tmp_path):
    """Loading a partial config should merge, not replace the entire config."""
    partial = {"parsing": {"strict_mode": True}}
    config_file = tmp_path / "partial.yaml"
    with open(config_file, "w") as f:
        yaml.dump(partial, f)

    cfg = Config()
    cfg.load_config(str(config_file))

    assert cfg.config["parsing"]["strict_mode"] is True
    assert "max_errors" in cfg.config["parsing"]
    assert "export" in cfg.config
    assert "translation" in cfg.config


def test_validate_config_negative_max_errors(capsys):
    cfg = Config()
    cfg.config["parsing"]["max_errors"] = -5
    cfg._validate_config()
    assert cfg.config["parsing"]["max_errors"] == 0
    captured = capsys.readouterr()
    assert "max_errors" in captured.out


def test_validate_config_negative_decimal_places(capsys):
    cfg = Config()
    cfg.config["output"]["decimal_places"] = -1
    cfg._validate_config()
    assert cfg.config["output"]["decimal_places"] == 2
    captured = capsys.readouterr()
    assert "decimal_places" in captured.out


def test_validate_config_negative_min_length(capsys):
    cfg = Config()
    cfg.config["orf"]["min_length"] = -10
    cfg._validate_config()
    assert cfg.config["orf"]["min_length"] == 0


def test_validate_config_negative_min_protein_length(capsys):
    cfg = Config()
    cfg.config["orf"]["min_protein_length"] = -3
    cfg._validate_config()
    assert cfg.config["orf"]["min_protein_length"] == 0


def test_validate_config_negative_default_k(capsys):
    cfg = Config()
    cfg.config["motif"]["default_k"] = -1
    cfg._validate_config()
    assert cfg.config["motif"]["default_k"] == 0


def test_validate_config_valid_values_unchanged():
    """Validation should not change valid values."""
    cfg = Config()
    cfg.config["parsing"]["max_errors"] = 50
    cfg.config["output"]["decimal_places"] = 4
    cfg.config["orf"]["min_length"] = 100
    cfg._validate_config()
    assert cfg.config["parsing"]["max_errors"] == 50
    assert cfg.config["output"]["decimal_places"] == 4
    assert cfg.config["orf"]["min_length"] == 100




def test_find_config_file_returns_none_when_no_file():
    cfg = Config()
    result = cfg._find_config_file()
    assert result is None or Path(result).exists()



def test_create_config_with_N_input(tmp_path, monkeypatch):
    """When user says N, should create config with defaults."""
    output_file = tmp_path / "config.yaml"
    cfg = Config()

    monkeypatch.setattr("builtins.input", lambda _: "N")
    cfg.create_config(str(output_file))

    assert output_file.exists()
    with open(output_file, "r") as f:
        saved = yaml.safe_load(f)
    assert saved["parsing"]["strict_mode"] is False
    assert saved["export"]["default_format"] == "csv"


def test_create_config_invalid_input_raises(tmp_path, monkeypatch):
    """When user enters invalid input, should raise ValueError."""
    output_file = tmp_path / "config.yaml"
    cfg = Config()

    monkeypatch.setattr("builtins.input", lambda _: "X")
    with pytest.raises(ValueError, match="Y or N"):
        cfg.create_config(str(output_file))


def test_create_config_overwrite_declined(tmp_path, monkeypatch, capsys):
    """When config exists and user declines overwrite, keep existing."""
    output_file = tmp_path / "config.yaml"
    output_file.write_text("existing: true")
    cfg = Config()

    monkeypatch.setattr("builtins.input", lambda _: "N")
    cfg.create_config(str(output_file))

    captured = capsys.readouterr()
    assert "keeping" in captured.out.lower()
    assert "existing" in output_file.read_text()

