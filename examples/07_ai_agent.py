"""M6 demo: the AI layer — natural language -> solved, costed flowsheet.

LLM backend is selectable and **local-first**. By default it uses a local model
via Ollama (no API key, no cost). Override with CALDYR_LLM_PROVIDER
(ollama | openai | anthropic) and CALDYR_LLM_MODEL.

    python examples/07_ai_agent.py                 # local Ollama (default)
    set CALDYR_LLM_MODEL=qwen3 & python examples/07_ai_agent.py
    set CALDYR_LLM_PROVIDER=openai & ...           # any OpenAI-compatible server

If no LLM backend is reachable, it runs the same tool calls a model would make
(a scripted transcript) so the tool layer is demonstrable offline.

Codex CLI users: register the MCP server instead and drive it from Codex:
    codex mcp add caldyr -- python -m caldyr.ai.mcp_server
"""
import os
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.ai import AgentSession, dispatch, ollama_available, run  # noqa: E402

REQUEST = ("Build me an ammonia synthesis loop and cost it. "
           "Report the levelized cost of ammonia per kg.")

# The tool calls a model makes (used as the offline demonstration).
AMMONIA_RECIPE = [
    ("new_flowsheet", {"components": ["nitrogen", "hydrogen", "ammonia", "argon"],
                       "property_package": "thermo:PR"}),
    ("add_unit", {"id": "MIX", "type": "Mixer", "params": {"dP": 0.0}}),
    ("add_unit", {"id": "PREHEAT", "type": "Heater", "params": {"T_out": 673.15}}),
    ("add_unit", {"id": "RXN", "type": "EquilibriumReactor",
                  "params": {"reaction": {"stoich": {"nitrogen": -1, "hydrogen": -3, "ammonia": 2},
                                          "key": "nitrogen"}, "T": 673.15}}),
    ("add_unit", {"id": "COOL", "type": "Heater", "params": {"T_out": 250.0}}),
    ("add_unit", {"id": "SEP", "type": "Flash", "params": {"T": 250.0, "P": 2e7}}),
    ("add_unit", {"id": "SPLIT", "type": "Splitter", "params": {"split": 0.9}}),
    ("add_feed", {"id": "MAKEUP", "to": "MIX:in1", "T": 300.0, "P": 2e7, "molar_flow": 100.0,
                  "z": {"nitrogen": 0.2475, "hydrogen": 0.7425, "argon": 0.01}}),
    ("connect", {"id": "S1", "from": "MIX:out", "to": "PREHEAT:in1"}),
    ("connect", {"id": "S2", "from": "PREHEAT:out", "to": "RXN:in1"}),
    ("connect", {"id": "S3", "from": "RXN:out", "to": "COOL:in1"}),
    ("connect", {"id": "S4", "from": "COOL:out", "to": "SEP:in1"}),
    ("connect", {"id": "PRODUCT", "from": "SEP:liquid", "to": None}),
    ("connect", {"id": "VAP", "from": "SEP:vapor", "to": "SPLIT:in1"}),
    ("connect", {"id": "RECYCLE", "from": "SPLIT:out1", "to": "MIX:in2"}),
    ("connect", {"id": "PURGE", "from": "SPLIT:out2", "to": None}),
    ("solve", {"tol": 1e-7}),
    ("cost", {"product_component": "ammonia"}),
]


def run_scripted() -> None:
    print("No LLM backend reachable — running the scripted tool transcript.\n")
    session = AgentSession()
    for name, args in AMMONIA_RECIPE:
        out = dispatch(session, name, args)
        flag = "ok " if out.get("ok") else "ERR"
        print(f"  [{flag}] {name:<14} {out.get('summary', out.get('error', ''))}")
    print("\nStart a local model (e.g. `ollama run qwen3`) to drive this from "
          "natural language instead.")


def main() -> None:
    provider = os.environ.get("CALDYR_LLM_PROVIDER", "ollama")
    print("=" * 64)
    print("  CALDYR AI LAYER - natural language -> solved, costed flowsheet")
    print("=" * 64)
    print(f'  provider: {provider}   request: "{REQUEST}"\n')

    if provider == "ollama" and not ollama_available():
        run_scripted()
        return

    res = run(REQUEST, provider=provider, verbose=True)
    print(f"\nDriven by {res.provider}:{res.model} - {len(res.tool_calls)} tool calls.")
    s = res.session
    if s.tea is not None:
        print(f"Result: LCOP ${s.tea.profitability.lcop:.3f}/kg ammonia, "
              f"TCI ${s.tea.capital.tci:,.0f}.")
    print("\n--- the model's answer ---")
    print(res.final_text)


if __name__ == "__main__":
    main()
