# AI copilot & MCP

Caldyr is AI-native: the engine exposes its actions as typed, JSON-schema'd
tools (`caldyr.ai.tools`), so a model operates the simulator exactly the way a
script does — it never invents numbers, the engine computes them.

## Tools

`new_flowsheet`, `add_unit`, `add_feed`, `connect`, `solve`, `cost`,
`optimize`, `stream_table`, `export_flow`, `list_unit_types`,
`list_property_packages`, `describe_flowsheet`, `explain_convergence`.

Errors return as data (`{ok: false, error}`), so models recover instead of
crashing.

## Providers — local-first

| Provider | Default model | Notes |
|---|---|---|
| `ollama` (default) | qwen3 | No key, no cost, fully local |
| `openai` | gpt-4o-mini | Any OpenAI-compatible server (LM Studio, vLLM...) |
| `anthropic` | claude-sonnet-4-6 | Opt-in only |

Select with `CALDYR_LLM_PROVIDER` / `CALDYR_LLM_MODEL`, or per call:
`caldyr.ai.run("build an ammonia loop and cost it", provider="ollama")`.

## In the web app

The Copilot panel talks to a persistent `ChatAgent` over `/ws/chat`: your
canvas flowsheet syncs into the session with each message, tool calls stream
back live, and any resulting flowsheet change is presented as a reviewable
diff. The one-click *Explain flowsheet* / *Diagnose solve* actions call the
diagnostic tools directly (`POST /ai/tool`) and need no LLM at all.

## MCP server

Use Caldyr from Claude Desktop, Codex CLI, or any MCP client:

```bash
codex mcp add caldyr -- python -m caldyr.ai.mcp_server
```

The same tool set is exposed over MCP, sharing one session per server process.
