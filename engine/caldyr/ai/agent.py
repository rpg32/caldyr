"""Natural-language agent over the Caldyr tools — local-first, provider-agnostic.

`run("build me an ammonia loop and cost it")` drives the typed engine tools in a
loop until the model is done. The LLM backend is selectable (default: a local
model via Ollama — no API key, no cost); see :mod:`caldyr.ai.llm`. The tool layer
itself is provider-independent and is what the tests exercise.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field

from .llm import LLMBackend, ToolCall, make_backend
from .session import AgentSession
from .tools import anthropic_tools, dispatch

SYSTEM_PROMPT = """You are Caldyr's flowsheet copilot. You build, solve, and cost \
chemical process flowsheets by calling the provided tools — you never invent \
numbers, the engine computes them.

EXACT PORT NAMES (use these literally as 'UNIT_ID:port'):
- Mixer: inlets in1, in2; outlet out.
- Heater (also used as a cooler): inlet in1; outlet out; energy duty.
- EquilibriumReactor / ConversionReactor: inlet in1; outlet out; energy duty.
- Flash: inlet in1; outlets vapor, liquid; energy duty.
- Splitter: inlet in1; outlets out1, out2.
- Pump / Compressor: inlet in1; outlet out; energy work.
- Valve: inlet in1; outlet out.

Rules:
- add_feed creates the feed stream AND wires it to the given port. NEVER call \
connect for a feed.
- connect is only for unit→unit streams and for products (set `to` to null). You \
do not need to connect energy/duty ports — leave them; the engine reports duties.
- Compositions are normalized for you; they need not sum to exactly 1.

Workflow:
1. new_flowsheet (thermo:PR for non-polar gases/hydrocarbons; thermo:NRTL for \
polar mixtures with azeotropes).
2. add_unit for each operation; add_feed for boundary feeds; connect unit→unit \
streams (and products to null).
3. solve, then cost (product_component = the sold species).
4. Briefly summarize: convergence, key flows, and economics (LCOP, capital, NPV).

Ammonia synthesis loop recipe: components [nitrogen, hydrogen, ammonia, argon], \
thermo:PR. Units: Mixer MIX; Heater PRE {T_out: 673.15}; EquilibriumReactor RXN \
{reaction: {stoich: {nitrogen: -1, hydrogen: -3, ammonia: 2}, key: nitrogen}, \
T: 673.15}; Heater COOL {T_out: 250}; Flash SEP {T: 250, P: 2e7}; Splitter SPLIT \
{split: 0.9}. Feed MAKEUP -> MIX:in1, 100 mol/s at 300 K and 2e7 Pa, \
z {nitrogen: 0.2475, hydrogen: 0.7425, argon: 0.01}. Connect MIX:out->PRE:in1, \
PRE:out->RXN:in1, RXN:out->COOL:in1, COOL:out->SEP:in1, SEP:vapor->SPLIT:in1, \
SPLIT:out1->MIX:in2 (recycle), SEP:liquid->null (product), SPLIT:out2->null \
(purge). Then solve and cost with product_component nitrogen->ammonia (use \
'ammonia'). Keep calling tools until solved AND costed, then stop and summarize."""


@dataclass
class AgentResult:
    final_text: str
    session: AgentSession
    tool_calls: list[str] = field(default_factory=list)
    turns: int = 0
    provider: str = ""
    model: str = ""


def run(prompt: str, *, provider: str | None = None, model: str | None = None,
        backend: LLMBackend | None = None, max_turns: int = 40,
        verbose: bool = False, **backend_opts) -> AgentResult:
    """Drive the tools from a natural-language request.

    provider: "ollama" (default, local), "openai" (OpenAI-compatible), or
    "anthropic". Pass a ready ``backend`` to override the factory entirely (e.g.
    a fake backend in tests).
    """
    llm = backend or make_backend(provider, **({"model": model} if model else {}), **backend_opts)
    session = AgentSession()
    tools = anthropic_tools()           # name/description/input_schema; backends adapt
    turns: list[dict] = [{"role": "user", "text": prompt}]
    called: list[str] = []
    final_text = ""

    for turn_i in range(1, max_turns + 1):
        resp = llm.complete(SYSTEM_PROMPT, turns, tools)
        if resp.text:
            final_text = resp.text
        turns.append({"role": "assistant", "text": resp.text, "tool_calls": resp.tool_calls})
        if not resp.tool_calls:
            break
        for call in resp.tool_calls:
            called.append(call.name)
            if verbose:
                print(f"  -> {call.name}({json.dumps(call.arguments)[:90]})")
            out = dispatch(session, call.name, call.arguments)
            if verbose:
                print(f"     {out.get('summary', out.get('error', ''))}")
            turns.append({"role": "tool", "tool_call_id": call.id, "name": call.name,
                          "content": json.dumps(out)})

    return AgentResult(final_text, session, called, turn_i,
                       provider=getattr(llm, "name", ""), model=getattr(llm, "model", ""))


__all__ = ["AgentResult", "run", "SYSTEM_PROMPT", "ToolCall"]
