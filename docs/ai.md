# AI copilot & MCP

Caldyr is AI-native: the engine exposes its actions as typed, JSON-schema'd
tools (`caldyr.ai.tools`), so a model operates the simulator exactly the way a
script does — it never invents numbers, the engine computes them.

## Tools

Build & edit: `new_flowsheet`, `add_unit`, `add_feed`, `connect`, `set_param`,
`remove_unit`, `remove_stream`.
Run: `solve`, `cost`, `optimize`.
Analyze: `stream_table`, `pinch_analysis` (heat-integration targets),
`sweep_parameter` (sensitivity studies on a copy — never alters the working
flowsheet), `property_table` (thermophysical properties over a T×P grid),
`describe_flowsheet`, `explain_convergence`.
Introspect: `list_unit_types` (pass `{type}` for one unit's full parameter
docs), `list_property_packages`, `export_flow`.

Errors return as data (`{ok: false, error}`), so models recover instead of
crashing.

## Providers — local-first

| Provider | Default model | Notes |
|---|---|---|
| `ollama` (default) | qwen3 | No key, no cost, fully local |
| `openai` | gpt-4o-mini | Any OpenAI-compatible server (LM Studio, vLLM...) |
| `anthropic` | claude-sonnet-4-6 | Opt-in only |

From Python, select with `CALDYR_LLM_PROVIDER` / `CALDYR_LLM_MODEL`, or per
call: `caldyr.ai.run("build an ammonia loop and cost it", provider="ollama")`.
Provider details:

- **Ollama** — install [Ollama](https://ollama.com), then `ollama pull qwen3`
  (or any tool-capable model). Point at a non-default host with `OLLAMA_HOST`.
- **OpenAI / compatible** — `OPENAI_API_KEY`, and `OPENAI_BASE_URL` to target a
  local server (LM Studio, vLLM, llama.cpp) instead of OpenAI.
- **Anthropic (Claude)** — install the extra (`pip install "caldyr[anthropic]"`)
  and set `ANTHROPIC_API_KEY` (or configure it in the app; see below).

## In the web app

The Copilot panel talks to a persistent `ChatAgent` over `/ws/chat`: your
canvas flowsheet syncs into the session with each message, tool calls stream
back live, and any resulting flowsheet change is presented as a reviewable
diff. The one-click *Explain flowsheet* / *Diagnose solve* actions call the
diagnostic tools directly (`POST /ai/tool`) and need no LLM at all.

### Configure your provider (bring your own LLM)

Open the Copilot panel and click the gear (or Ctrl+K → "Copilot: AI provider
settings"). Pick **Ollama** (local, no key), **OpenAI / compatible** (your key
or a local server's base URL), or **Anthropic** (your Claude key), then **Test
connection** and **Save**. The setting — including any API key — is stored by
the local Caldyr server (`~/.config/caldyr/llm.json`, or `%APPDATA%\caldyr` on
Windows), **not in the browser**; the key is sent once over the loopback API and
never kept client-side. The panel shows which providers are reachable and warns
before you send a message if none is ready.

## MCP server

Use Caldyr from Claude Desktop, Codex CLI, or any MCP client:

```bash
codex mcp add caldyr -- python -m caldyr.ai.mcp_server
```

The same tool set is exposed over MCP, sharing one session per server process.
