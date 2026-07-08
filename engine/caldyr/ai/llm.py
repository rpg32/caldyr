"""Pluggable LLM backends for the agent — **local first**.

Development defaults to a local model via **Ollama** (no API key, no cost, no data
leaving the machine). An OpenAI-compatible backend covers other local servers
(LM Studio, vLLM, llama.cpp, ...) and hosted OpenAI; an Anthropic backend is
available but opt-in. Codex CLI users connect via the MCP server
(:mod:`caldyr.ai.mcp_server`) instead of a chat backend.

The agent speaks a neutral conversation format (lists of role dicts + ToolCall
objects); each backend translates it to its own wire protocol, so the agent loop
never depends on a provider's message schema.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Protocol, cast

import httpx

def _normalize_host(host: str) -> str:
    """Make a connectable client URL: add a scheme if missing, and map the
    server bind address 0.0.0.0 to loopback (OLLAMA_HOST is often a bind addr)."""
    host = (host or "").strip()
    if not host.startswith(("http://", "https://")):
        host = "http://" + host
    return host.replace("0.0.0.0", "127.0.0.1").rstrip("/")


DEFAULT_OLLAMA_HOST = _normalize_host(os.environ.get("OLLAMA_HOST") or "http://localhost:11434")
_TIMEOUT = httpx.Timeout(600.0)


@dataclass
class ToolCall:
    id: str
    name: str
    arguments: dict


@dataclass
class LLMResponse:
    text: str = ""
    tool_calls: list[ToolCall] = field(default_factory=list)


# A conversation turn is one of:
#   {"role": "user", "text": str}
#   {"role": "assistant", "text": str, "tool_calls": list[ToolCall]}
#   {"role": "tool", "tool_call_id": str, "name": str, "content": str}
Turn = dict


class LLMBackend(Protocol):
    name: str
    model: str
    def complete(self, system: str, turns: list[Turn], tools: list[dict]) -> LLMResponse: ...


def _openai_tool_schema(tools: list[dict]) -> list[dict]:
    """Caldyr tool schema -> OpenAI/Ollama function-tool schema."""
    return [{"type": "function",
             "function": {"name": t["name"], "description": t["description"],
                          "parameters": t["input_schema"]}}
            for t in tools]


# -- Ollama (default, local) -----------------------------------------------
class OllamaBackend:
    """Local models via Ollama's native /api/chat tool-calling."""
    name = "ollama"

    def __init__(self, model: str | None = None, host: str = DEFAULT_OLLAMA_HOST,
                 think: bool = False) -> None:
        self.model = model or os.environ.get("CALDYR_LLM_MODEL") or "qwen3"
        self.host = _normalize_host(host)
        self.think = think

    def complete(self, system: str, turns: list[Turn], tools: list[dict]) -> LLMResponse:
        messages = [{"role": "system", "content": system}]
        for t in turns:
            if t["role"] == "user":
                messages.append({"role": "user", "content": t["text"]})
            elif t["role"] == "assistant":
                m: dict = {"role": "assistant", "content": t.get("text", "")}
                if t.get("tool_calls"):
                    m["tool_calls"] = [{"id": c.id, "type": "function",
                                        "function": {"name": c.name, "arguments": c.arguments}}
                                       for c in t["tool_calls"]]
                messages.append(m)
            else:  # tool result
                messages.append({"role": "tool", "content": t["content"],
                                 "tool_name": t["name"]})
        body = {"model": self.model, "messages": messages,
                "tools": _openai_tool_schema(tools), "stream": False, "think": self.think}
        data = httpx.post(f"{self.host}/api/chat", json=body, timeout=_TIMEOUT).json()
        if "error" in data:
            raise RuntimeError(f"ollama error: {data['error']}")
        msg = data.get("message", {})
        calls = []
        for i, tc in enumerate(msg.get("tool_calls", []) or []):
            fn = tc.get("function", {})
            args = fn.get("arguments", {})
            if isinstance(args, str):
                args = json.loads(args or "{}")
            calls.append(ToolCall(id=tc.get("id") or f"call_{i}", name=fn["name"], arguments=args))
        return LLMResponse(text=msg.get("content", "") or "", tool_calls=calls)


# -- OpenAI-compatible (other local servers or hosted OpenAI) ---------------
class OpenAIBackend:
    name = "openai"

    def __init__(self, model: str | None = None, base_url: str | None = None,
                 api_key: str | None = None) -> None:
        self.model = model or os.environ.get("CALDYR_LLM_MODEL") or "gpt-4o-mini"
        self.base_url = (base_url or os.environ.get("OPENAI_BASE_URL")
                         or "https://api.openai.com/v1").rstrip("/")
        self.api_key = api_key or os.environ.get("OPENAI_API_KEY") or "not-needed"

    def complete(self, system: str, turns: list[Turn], tools: list[dict]) -> LLMResponse:
        messages = [{"role": "system", "content": system}]
        for t in turns:
            if t["role"] == "user":
                messages.append({"role": "user", "content": t["text"]})
            elif t["role"] == "assistant":
                m: dict = {"role": "assistant", "content": t.get("text", "") or None}
                if t.get("tool_calls"):
                    m["tool_calls"] = [
                        {"id": c.id, "type": "function",
                         "function": {"name": c.name, "arguments": json.dumps(c.arguments)}}
                        for c in t["tool_calls"]]
                messages.append(m)
            else:
                messages.append({"role": "tool", "tool_call_id": t["tool_call_id"],
                                 "content": t["content"]})
        body = {"model": self.model, "messages": messages,
                "tools": _openai_tool_schema(tools)}
        headers = {"Authorization": f"Bearer {self.api_key}"}
        data = httpx.post(f"{self.base_url}/chat/completions", json=body,
                          headers=headers, timeout=_TIMEOUT).json()
        if "error" in data:
            raise RuntimeError(f"openai backend error: {data['error']}")
        msg = data["choices"][0]["message"]
        calls = []
        for tc in msg.get("tool_calls") or []:
            fn = tc["function"]
            calls.append(ToolCall(id=tc["id"], name=fn["name"],
                                  arguments=json.loads(fn.get("arguments") or "{}")))
        return LLMResponse(text=msg.get("content") or "", tool_calls=calls)


# -- Anthropic (opt-in) -----------------------------------------------------
class AnthropicBackend:
    name = "anthropic"

    def __init__(self, model: str | None = None, max_tokens: int = 4096,
                 api_key: str | None = None, base_url: str | None = None) -> None:
        self.model = model or os.environ.get("CALDYR_LLM_MODEL") or "claude-sonnet-4-6"
        self.max_tokens = max_tokens
        # Fall through to the SDK's own ANTHROPIC_API_KEY env read when unset.
        self.api_key = api_key or os.environ.get("ANTHROPIC_API_KEY") or None
        self.base_url = base_url or os.environ.get("ANTHROPIC_BASE_URL") or None

    def complete(self, system: str, turns: list[Turn], tools: list[dict]) -> LLMResponse:
        import anthropic

        client_opts: dict = {}
        if self.api_key:
            client_opts["api_key"] = self.api_key
        if self.base_url:
            client_opts["base_url"] = self.base_url
        messages = []
        for t in turns:
            if t["role"] == "user":
                messages.append({"role": "user", "content": t["text"]})
            elif t["role"] == "assistant":
                content: list[dict] = []
                if t.get("text"):
                    content.append({"type": "text", "text": t["text"]})
                for c in t.get("tool_calls", []):
                    content.append({"type": "tool_use", "id": c.id, "name": c.name,
                                    "input": c.arguments})
                messages.append({"role": "assistant", "content": content})
            else:
                messages.append({"role": "user", "content": [
                    {"type": "tool_result", "tool_use_id": t["tool_call_id"],
                     "content": t["content"]}]})
        atools = [{"name": t["name"], "description": t["description"],
                   "input_schema": t["input_schema"]} for t in tools]
        resp = anthropic.Anthropic(**client_opts).messages.create(
            model=self.model, system=system, max_tokens=self.max_tokens,
            tools=atools, messages=messages)  # type: ignore[arg-type]
        text = "".join(b.text for b in resp.content if b.type == "text")
        calls = [ToolCall(id=b.id, name=b.name, arguments=cast("dict", b.input))
                 for b in resp.content if b.type == "tool_use"]
        return LLMResponse(text=text, tool_calls=calls)


_BACKENDS = {"ollama": OllamaBackend, "openai": OpenAIBackend, "anthropic": AnthropicBackend}


def make_backend(provider: str | None = None, **opts) -> LLMBackend:
    """Build the selected backend. Defaults to ``CALDYR_LLM_PROVIDER`` or
    ``"ollama"`` — local-first, no key, no cost."""
    provider = (provider or os.environ.get("CALDYR_LLM_PROVIDER") or "ollama").lower()
    if provider not in _BACKENDS:
        raise ValueError(f"unknown LLM provider {provider!r}; "
                         f"choose from {sorted(_BACKENDS)}")
    return _BACKENDS[provider](**opts)


def ollama_available(host: str = DEFAULT_OLLAMA_HOST) -> bool:
    try:
        return httpx.get(f"{_normalize_host(host)}/api/tags", timeout=2.0).status_code == 200
    except Exception:
        return False


def anthropic_available() -> bool:
    """Whether the optional ``anthropic`` SDK is importable (the ``.[anthropic]``
    extra); the Anthropic backend needs it."""
    try:
        import anthropic  # noqa: F401
        return True
    except Exception:
        return False


def ollama_tool_models(host: str = DEFAULT_OLLAMA_HOST) -> list[str]:
    """Names of locally-pulled Ollama models that support tool calling."""
    try:
        tags = httpx.get(f"{_normalize_host(host)}/api/tags", timeout=3.0).json()
    except Exception:
        return []
    out = []
    for m in tags.get("models", []):
        caps = m.get("details", {}).get("families") or []
        if "tools" in (m.get("capabilities") or []) or caps:
            out.append(m["name"])
    return out
