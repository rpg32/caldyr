# Getting started

## Install

```bash
git clone <repo> caldyr && cd caldyr
pip install -e ".[dev,api,ai,kinetics]"   # engine + API + AI + Cantera
pytest -q                                  # the full validated suite
```

The extras: `dev` (pytest/ruff/mypy), `api` (FastAPI bridge), `ai` (local-LLM
copilot + MCP server), `kinetics` (Cantera, for the Gibbs reactor). Also
available: `eo` (Pyomo/IDAES equation-oriented backend) and `anthropic`
(hosted LLM backend, opt-in only).

## Run the app

```bash
# 1. engine API
python -m uvicorn api.main:app --port 8753

# 2. web UI (separate terminal)
cd web && npm install && npm run dev      # http://localhost:5273
```

Optional, for the AI copilot: install [Ollama](https://ollama.com) and pull a
tool-calling model (`ollama pull qwen3`). No API keys needed.

## First flowsheet (60 seconds)

1. Open **Projects → Templates → Ammonia loop** (or one of the book plants —
   cyclohexane, VCM, DME; or build from scratch: add components in the
   inspector, drag units from the palette, connect ports).
2. Press **Solve** — watch live convergence in the status bar; the Streams tab
   shows the table, convergence plot, and balance check.
3. Pick a product component, press **Cost** — LCOP, capital, opex, NPV, tornado
   sensitivity, and optional Monte-Carlo.
4. Open the **Copilot** and ask it to change something — its edit arrives as a
   diff you accept or reject.

## Headless (Python)

```python
from caldyr.io import load_flow
fs = load_flow("flowsheet.flow")
report = fs.solve()                        # or backend="equation_oriented"

from caldyr.economics import TEAConfig, analyze
tea = analyze(fs, report, TEAConfig(product_component="ammonia"))
print(tea.profitability.lcop)
```

Every example in `examples/` (01–18) is runnable and validated — the
authoritative feature tour, from a two-unit balance through rigorous columns,
absorbers, extraction, pipe networks & relief sizing, heat integration, and
the full cyclohexane plant of the book's ch. 15. `examples/04`–`05` build and
cost the ammonia loop end-to-end.

## This documentation site

```bash
python -m pip install --user mkdocs-material
python -m mkdocs serve        # from the repo root
```
