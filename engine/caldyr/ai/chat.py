"""Multi-turn chat over the Caldyr tools — the engine side of the web chat panel.

Unlike :func:`caldyr.ai.agent.run` (one-shot), a :class:`ChatAgent` keeps the
conversation and the working flowsheet across messages, can be synced with an
external flowsheet (the canvas) before each message, and reports every step
through an event callback so a UI can stream progress:

    {"type": "text",        "text": str}                       # assistant prose
    {"type": "tool_call",   "name": str, "arguments": dict}
    {"type": "tool_result", "name": str, "ok": bool, "summary": str}

The final return carries the resulting flowsheet as a ``.flow`` dict (or None)
so the UI can diff it against the canvas and offer accept/reject.
"""
from __future__ import annotations

import json
from typing import Any, Callable

from ..io import from_dict, to_dict
from .agent import SYSTEM_PROMPT
from .llm import LLMBackend, make_backend
from .session import AgentSession
from .tools import anthropic_tools, dispatch

OnEvent = Callable[[dict], None]

CHAT_SYSTEM = SYSTEM_PROMPT + """

You are running inside the Caldyr web app as a chat copilot. The user's canvas
flowsheet is loaded into your session before each message — edit THAT flowsheet
(do not call new_flowsheet unless the user asks to start over). After structural
edits, solve before answering questions about stream values. Use
describe_flowsheet / explain_convergence to ground explanations. Keep final
answers short and concrete; the UI shows the numbers, you provide the insight."""


class ChatAgent:
    """A persistent conversation bound to one AgentSession."""

    def __init__(self, *, provider: str | None = None, model: str | None = None,
                 backend: LLMBackend | None = None, **backend_opts: Any) -> None:
        self.llm: LLMBackend = backend or make_backend(
            provider, **({"model": model} if model else {}), **backend_opts)
        self.session = AgentSession()
        self.turns: list[dict] = []
        self.tools = anthropic_tools()

    # -- canvas sync --------------------------------------------------------
    def load_flow(self, doc: dict) -> None:
        """Replace the session flowsheet with the canvas state (pre-message sync)."""
        self.session.flowsheet = from_dict(doc)
        self.session.report = None
        self.session.tea = None

    def export_flow(self) -> dict | None:
        if self.session.flowsheet is None:
            return None
        return to_dict(self.session.flowsheet)

    # -- one user message ----------------------------------------------------
    def send(self, text: str, *, on_event: OnEvent | None = None,
             max_steps: int = 30) -> dict:
        """Process one user message; returns {text, flow, tool_calls}."""
        emit = on_event or (lambda _e: None)
        self.turns.append({"role": "user", "text": text})
        called: list[str] = []
        final_text = ""

        for _step in range(max_steps):
            resp = self.llm.complete(CHAT_SYSTEM, self.turns, self.tools)
            if resp.text:
                final_text = resp.text
                emit({"type": "text", "text": resp.text})
            self.turns.append({"role": "assistant", "text": resp.text,
                               "tool_calls": resp.tool_calls})
            if not resp.tool_calls:
                break
            for call in resp.tool_calls:
                called.append(call.name)
                emit({"type": "tool_call", "name": call.name, "arguments": call.arguments})
                out = dispatch(self.session, call.name, call.arguments)
                emit({"type": "tool_result", "name": call.name,
                      "ok": bool(out.get("ok", True)),
                      "summary": str(out.get("summary") or out.get("error") or "")[:400]})
                self.turns.append({"role": "tool", "tool_call_id": call.id,
                                   "name": call.name, "content": json.dumps(out)})

        return {"text": final_text, "flow": self.export_flow(), "tool_calls": called}


__all__ = ["ChatAgent", "CHAT_SYSTEM"]
