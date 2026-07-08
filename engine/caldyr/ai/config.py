"""Persisted LLM configuration for the web copilot.

An installed user picks their provider/model/key once in the app; it is stored
**server-side** (this file) rather than in the browser, so the API key never
lives in localStorage (readable by any page script) or in the process
environment. The server is loopback-only and single-user, so a 0600 file in the
user config directory is an appropriate store; swapping in the OS keyring later
is a drop-in change behind ``load_config`` / ``save_config``.

Resolution order for any field is: saved config -> environment variable ->
backend default (the backends in :mod:`caldyr.ai.llm` already apply the last
two, so unset fields simply fall through).
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

_FIELDS = ("provider", "model", "base_url", "api_key")


def config_dir() -> Path:
    """Per-user config directory (``%APPDATA%\\caldyr`` on Windows, else
    ``$XDG_CONFIG_HOME/caldyr`` or ``~/.config/caldyr``)."""
    if os.name == "nt":
        base = os.environ.get("APPDATA") or str(Path.home() / "AppData" / "Roaming")
    else:
        base = os.environ.get("XDG_CONFIG_HOME") or str(Path.home() / ".config")
    return Path(base) / "caldyr"


def config_path() -> Path:
    return config_dir() / "llm.json"


def load_config() -> dict[str, Any]:
    """Return the saved config, or ``{}`` if none/unreadable."""
    try:
        data = json.loads(config_path().read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return {k: v for k, v in data.items() if k in _FIELDS and v not in (None, "")}


def save_config(cfg: dict[str, Any]) -> dict[str, Any]:
    """Persist provider/model/base_url/api_key (only the known, non-empty
    fields). Written atomically with owner-only permissions. Returns the
    saved config (sans key) for convenience."""
    clean = {k: cfg[k] for k in _FIELDS if cfg.get(k) not in (None, "")}
    d = config_dir()
    d.mkdir(parents=True, exist_ok=True)
    path = config_path()
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(clean, indent=2), encoding="utf-8")
    try:  # POSIX: owner read/write only. No-op-ish on Windows, harmless.
        os.chmod(tmp, 0o600)
    except OSError:
        pass
    os.replace(tmp, path)
    return public_config(clean)


def public_config(cfg: dict[str, Any] | None = None) -> dict[str, Any]:
    """The config safe to hand back to the browser: never includes the key,
    only whether one is set."""
    c = dict(load_config() if cfg is None else cfg)
    return {
        "provider": c.get("provider"),
        "model": c.get("model"),
        "base_url": c.get("base_url"),
        "has_key": bool(c.get("api_key")),
    }


def backend_kwargs(overrides: dict[str, Any] | None = None) -> tuple[str | None, dict[str, Any]]:
    """Resolve saved config (with optional per-call overrides) into
    ``(provider, opts)`` suitable for :func:`caldyr.ai.llm.make_backend`.

    Maps ``base_url`` to the per-provider argument name (``host`` for Ollama,
    ``base_url`` otherwise) and drops empty values so the backend's own env /
    default resolution still applies to anything the user did not set."""
    cfg = {**load_config(), **{k: v for k, v in (overrides or {}).items() if v not in (None, "")}}
    provider = (cfg.get("provider") or None)
    prov = (provider or "").lower()
    opts: dict[str, Any] = {}
    if cfg.get("model"):
        opts["model"] = cfg["model"]
    if cfg.get("base_url"):
        opts["host" if prov == "ollama" else "base_url"] = cfg["base_url"]
    if cfg.get("api_key") and prov != "ollama":
        opts["api_key"] = cfg["api_key"]
    return provider, opts
