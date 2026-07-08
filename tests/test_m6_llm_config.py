"""BYO-LLM server-side config store (caldyr.ai.config)."""
import importlib

import pytest


@pytest.fixture
def cfg(tmp_path, monkeypatch):
    """Point the config dir at a temp location so tests never touch the real
    user config, and reload the module so config_dir() re-reads the env."""
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    monkeypatch.setenv("APPDATA", str(tmp_path))  # Windows path
    import caldyr.ai.config as c
    importlib.reload(c)
    return c


def test_roundtrip_and_key_is_hidden(cfg):
    assert cfg.load_config() == {}
    pub = cfg.save_config({"provider": "anthropic", "model": "claude-x",
                           "api_key": "sk-secret"})
    # public view never leaks the key
    assert pub == {"provider": "anthropic", "model": "claude-x",
                   "base_url": None, "has_key": True}
    # raw load has it (server-side use only)
    assert cfg.load_config()["api_key"] == "sk-secret"


def test_empty_fields_are_dropped(cfg):
    cfg.save_config({"provider": "ollama", "model": "", "base_url": None,
                     "api_key": ""})
    assert cfg.load_config() == {"provider": "ollama"}


def test_backend_kwargs_maps_base_url_per_provider(cfg):
    cfg.save_config({"provider": "ollama", "base_url": "http://host:1234"})
    provider, opts = cfg.backend_kwargs()
    assert provider == "ollama"
    assert opts == {"host": "http://host:1234"}          # ollama uses host

    cfg.save_config({"provider": "openai", "base_url": "http://x/v1",
                     "api_key": "k", "model": "m"})
    provider, opts = cfg.backend_kwargs()
    assert provider == "openai"
    assert opts == {"model": "m", "base_url": "http://x/v1", "api_key": "k"}


def test_backend_kwargs_never_sends_key_to_ollama(cfg):
    cfg.save_config({"provider": "ollama", "api_key": "leak"})
    _, opts = cfg.backend_kwargs()
    assert "api_key" not in opts


def test_overrides_take_precedence(cfg):
    cfg.save_config({"provider": "ollama", "model": "qwen3"})
    provider, opts = cfg.backend_kwargs({"model": "other"})
    assert provider == "ollama"
    assert opts["model"] == "other"
