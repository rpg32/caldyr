# Caldyr

**Open-source, AI-native process simulation & techno-economic analysis.**
Build BFDs / PFDs / P&IDs, run steady-state mass & energy balances, size and
cost equipment, and optimize and scale flowsheets — from a modern web UI *and* a
scriptable Python API.

> Working name; rename freely. Think *open Aspen HYSYS + economic analyzer, but
> git-native, scriptable, and with an LLM in the loop.*

## Why
Commercial simulators (Aspen HYSYS/Plus, PRO/II, UniSim, gPROMS; Thermoflow for
power) are costly, Windows-locked, GUI-only, black-box, and unfriendly to version
control, automation, and optimization. Caldyr is **open, free, API-first,
git-native, AI-native**, with **techno-economics and optimization in the core** —
and it stands on validated open thermodynamics rather than reinventing it.

## What works today
- **18 unit operations**, each validated against a cited reference: mixers,
  heaters/coolers, fired equipment, splitters, valves, pumps, compressors,
  expanders, flash drums, three-phase (VLLE) separators, component splitters,
  heat exchangers, conversion / equilibrium / Gibbs (Cantera) / CSTR / PFR
  reactors, and both shortcut (FUG) and rigorous (MESH) distillation columns.
- **Two solver backends, one physics** — sequential-modular (Wegstein tearing)
  and equation-oriented — agreeing to ~1e-9; optimization on top.
- **Logical ops** (Set/Adjust), solver hints, and balance diagnostics.
- **Techno-economics in core**: Turton costing, CEPCI escalation, capital/opex,
  NPV/IRR/LCOP, tornado sensitivity, Monte-Carlo uncertainty.
- **Web app**: flowsheet canvas with BFD/PFD/P&ID views, phase/temperature
  stream coloring, auto-layout, groups, undo/copy/paste, projects & templates,
  optimization & case-study builders, plots — and an **AI copilot** (local LLM
  via Ollama) whose every edit arrives as a reviewable diff.
- **MCP server** so Claude Desktop / Codex CLI can drive the engine directly.

## Foundations (we build on, not reinvent)
`thermo` + `chemicals` + `CoolProp` (properties/EOS) · `Cantera` (equilibrium/
kinetics) · `Pyomo` (optional EO optimization backend) · `@xyflow/react`
(canvas) · DWSIM (interop/validation oracle, GPL — API only).

## Quickstart
```bash
pip install -e ".[dev,api,ai,kinetics]"
pytest -q                                  # the full validated suite

python -m uvicorn api.main:app --port 8753 # engine API
cd web && npm install && npm run dev       # UI at http://localhost:5273
```
Docs: `python -m mkdocs serve` (or read `docs/`). Start with
`docs/getting-started.md`; architecture in `ARCHITECTURE.md`.

## Layout
```
engine/caldyr/{core,thermo,unitops,solver,economics,analysis,ai,io}   # the engine
api/   FastAPI bridge        web/  React + @xyflow/react UI
examples/  validated cases   tests/  pytest    docs/  mkdocs site
```

## License
Permissive (see `LICENSE`). Never vendor GPL code (e.g. DWSIM internals) — interop
across its API only.
