# CLAUDE.md — Caldyr

> **Caldyr** (working name — rename freely) is an open-source, AI-native process
> simulation and techno-economic analysis platform: build BFDs / PFDs / P&IDs,
> run steady-state mass & energy balances, size and cost equipment, and optimize
> and scale flowsheets — all from a modern web UI *and* a scriptable Python API.
> Think "open Aspen HYSYS + Aspen Process Economic Analyzer, but git-native,
> scriptable, and with an LLM in the loop."

This file orients you (Claude Code) to the project. Read `ARCHITECTURE.md`,
`ROADMAP.md`, and `docs/DATA_MODEL.md` before writing code.

---

## 1. Mission & wedge

Commercial process simulators (Aspen HYSYS/Plus, PRO/II, UniSim, ProMax,
gPROMS; Thermoflow/EBSILON for power) are expensive, Windows-locked, GUI-only,
black-box, and hostile to version control and automation. They are also weak at
*equation-oriented optimization* and *integrated techno-economics*.

**We do not compete on thermodynamics.** Validated property methods already
exist in open form. We compete on:

1. **Open & free** (permissive license) — no seat licenses.
2. **Git-native** — flowsheets are plain JSON (`.flow`), diffable and reviewable.
3. **API-first** — every action available in Python; the GUI is a thin layer
   over the engine. Flowsheet-as-code and flowsheet-as-canvas are equivalent.
4. **AI-native** — Claude can build, edit, explain, debug-convergence, and
   report on flowsheets through a typed tool interface over the engine.
5. **Techno-economics as a first-class citizen** — capex/opex/NPV/scaling live
   in the core, not a bolted-on separate product.
6. **Optimization & scale-up built in** — via an IDAES/Pyomo backend, the thing
   Aspen does worst.

## 2. Build on proven foundations (do NOT reinvent these)

| Need | Use | Notes |
|------|-----|-------|
| Pure & mixture properties, EOS (PR, SRK, etc.) | `thermo` + `chemicals` (Caleb Bell), `CoolProp` | MIT-licensed, validated |
| Equilibrium / kinetics | `Cantera` | BSD |
| Equation-oriented modeling, optimization, scale-up | `IDAES` + `Pyomo` (+ IPOPT) | DOE/NETL, BSD-style; our heavyweight backend |
| Numerics | `numpy`, `scipy` | |
| Flowsheet canvas | `@xyflow/react` (React Flow v12) | TS/React |
| Validation oracle | DWSIM (headless TCP/IP server + Python API) | GPLv3 — **interop only, never vendor its code** |

**License hygiene:** Caldyr is permissively licensed (see `LICENSE`). Never copy
or statically link GPL code (e.g., DWSIM internals). Interop with DWSIM only
across its API/process boundary.

## 3. Architecture in one breath

A layered engine (`engine/caldyr/`):
`core` (data model: Component, Stream, UnitOp, Flowsheet) →
`thermo` (PropertyPackage interface wrapping thermo/CoolProp/Cantera) →
`unitops` (UnitOp models implementing a common solve contract) →
`solver` (sequential-modular first: topological sort, tear-stream detection,
recycle convergence; IDAES/Pyomo equation-oriented backend later) →
`economics` (sizing → costing → profitability → scaling) →
`io` (`.flow` JSON load/save).
`api/` exposes the engine over FastAPI; `web/` is the React/`@xyflow/react` UI;
the AI layer exposes engine actions as typed tools. Full detail in `ARCHITECTURE.md`.

## 4. Core data model (see docs/DATA_MODEL.md)

- **Component** — a chemical species + property metadata.
- **Stream** — material stream: T, P, molar flow, composition {component: mole frac},
  resolved phase/enthalpy from the property package. Plus energy streams.
- **Port** — a typed connection point on a unit op (inlet/outlet, material/energy).
- **UnitOp** — abstract; holds parameters, ports; implements `solve(inlets) -> outlets`.
- **Flowsheet** — directed graph of UnitOps connected by Streams; owns the
  component list and the selected PropertyPackage.

## 5. Conventions

**Python** (engine/api):
- 3.11+, full type hints, `dataclasses`/`pydantic` for models.
- Format/lint with `ruff`; type-check with `mypy` (or `pyright`).
- Tests with `pytest`. **Every unit op and property method ships with a test
  validated against a textbook or DWSIM/Aspen reference case** (cite the source
  in the test).
- SI units internally (K, Pa, mol/s, J). Conversions only at the I/O edges.
- No silent failures in solvers: raise typed exceptions with diagnostics.

**TypeScript** (web): strict mode, functional React components, `@xyflow/react`
for the canvas. The web app never contains physics — it calls the engine API.

**General:** small, reviewable commits. Conventional Commits. Keep the engine
importable and runnable headless without the web app.

## 6. Build / test commands

```bash
# engine
cd engine && uv pip install -e ".[dev]"   # or: pip install -e ".[dev]"
pytest -q
ruff check . && mypy caldyr

# api
uvicorn api.main:app --reload

# web (once scaffolded)
cd web && npm install && npm run dev
```

## 7. Current milestone — M7 (see ROADMAP.md for the full plan)

**M0–M6 ✅ are complete and green** (66 passing tests; ruff + mypy clean; web app
builds + verified in-browser).
What exists today:
- `core` (Component, Stream/EnergyStream, Port, UnitOp, Flowsheet), `io` (`.flow`
  exact round-trip).
- `thermo`: two backends behind the `PropertyPackage` protocol, selected by the
  flowsheet's `property_package` string — `ThermoPackage` (`thermo:PR`/`thermo:SRK`,
  cubic EOS, for non-polar systems) and `ActivityPackage` (`thermo:NRTL`, gamma-phi
  with ChemSep parameters, for polar mixtures/azeotropes). Both expose enthalpy,
  entropy, volume, PT/PH/PS flash, bubble/dew, per-phase VLE, and `lnKeq` (reaction
  equilibrium constant). **Enthalpy is formation-inclusive** (absolute basis), so
  reactor heats of reaction and adiabatic temperatures are correct, while
  conserved-composition balances still close (the offset cancels). Phase
  identification is by molar volume (robust at high P where thermo mis-labels
  vapor vs liquid). NRTL validated vs the ethanol/water azeotrope; `lnKeq` vs the
  Haber-Bosch ammonia equilibrium.
- `unitops`: Mixer, Heater, Splitter, Valve, Pump, Compressor, FlashDrum,
  HeatExchanger (LMTD + effectiveness-NTU), MultiStreamExchanger (N-pass LNG/
  plate-fin; weighted zone composite-curve analysis, MITA/UA, global
  min_approach/ua specs — Hameed §9.5.2), Saturator (gas humidifier / HYSYS
  Stream Saturator analogue — saturation y* from a probe flash, RH knob, latent
  duty), ConversionReactor, EquilibriumReactor,
  ShortcutColumn (FUG: Fenske/Underwood/Gilliland/Kirkbride, validated vs Wankat;
  sized & costed as tower + sieve trays + condenser/reboiler), Expander (isentropic
  turbine, Turton axial-turbine costing), ComponentSplitter (black-box splits with
  honest duty), ThreePhaseSeparator (VLLE via thermo FlashVLN `flash_pt_3p`; PR/SRK
  only; requires a T spec — PH 3-phase deliberately unsupported), GibbsReactor
  (Cantera equilibrate over mapped flowsheet species only), CSTR + PFR (power-law
  kinetics, validated vs Fogler/Levenspiel closed forms), RigorousColumn
  (Wang–Henke bubble-point MESH, full stage profiles on `unit.design`,
  validated vs the FUG design point; absorber mode not yet), FiredHeater
  (fuel duty = Q/η, fired-fuel opex), AirCooler (approach-validated, fan
  electricity opex).
- `analysis/pinch.py`: problem-table pinch targeting (validated vs Kemp 2e),
  composite curves, recovery potential — the heat-integration tool.
- `solver`: third backend `"pyomo"` (pynumero grey-box over the same residuals;
  live solve needs cyipopt — no Windows wheels, conda flips the skip).
- economics sizers are a registry now (`@register_sizer`); versions at 0.9.0;
  CI in `.github/workflows/ci.yml`; suite ~233 tests + stress flowsheet.
- `core/components_db.py`: `resolve_component` via `chemicals` + 120-species
  curated catalog (served at `GET /components`); economics MW lookups fall back
  to chemicals (no more KeyError on uncommon species).
- `caldyr.ai` additions: `diagnostics.py` (describe_flowsheet,
  explain_convergence — registered tools), `chat.py` (`ChatAgent`, multi-turn,
  canvas-synced; drives the web Copilot over the API's `/ws/chat`).
- API (Phase C/D additions): `/envelope`, `monte_carlo` on `/cost`,
  `report.history`, `/balance`, `/ai/tool`, `/components`, WS `/ws/solve`
  (live residuals via the solver's `on_iteration` callback) and `/ws/chat`.
- logical ops (flowsheet-level, in `.flow` under `"logical"`): Set (param binding)
  and Adjust (spec-driven outer root-find, backend-agnostic), plus `solver_hints`
  (tear guesses/tolerance) and `solver.balance_report(fs)` diagnostics.
  Reactions are JSON-friendly dicts (`{"stoich": {...}, "key": ...}`) so `.flow`
  round-trips. Each unit validated against first principles / textbook data.
- `solver`: two backends via `fs.solve(backend=...)`. `"sequential"` —
  `SequentialModularSolver`, acyclic single-sweep fast path; recycles torn via DFS
  back-edges and converged by Wegstein. `"equation_oriented"` —
  `EquationOrientedSolver`, assembles all unit equations into one scaled residual
  and solves simultaneously with scipy Newton (recycles need no tearing); reuses
  the same `unit.solve`/property physics, so the two backends agree to ~1e-9.
  `solver.optimize(fs, objective, design_vars, constraints)` minimizes any
  callable of the solved flowsheet over unit parameters (scipy SLSQP, solve in the
  loop, auto-scaled objective). `SolveReport` carries tear streams, residual
  history, tolerance, method.
- Examples: `01_mixer_heater` (M0), `02_flash_recycle` (M1), `03_azeotrope_nrtl`
  (NRTL vs PR), `04_ammonia_loop` (M2 DoD — a converging Haber-Bosch loop with
  reactor + separator + recycle + inert purge; atom balances close).

- `economics`: full TEA pipeline behind `economics.analyze(fs, report, TEAConfig)`
  → sizing (`sizing.py`, area/power/volume + utility selection) → bare-module
  costing (`costing.py`, Turton 4e correlations in typed `data.py` tables, CEPCI
  escalation, parallel trains, vessel wall-stress Fp) → capital (`capital.py`,
  ISBL/OSBL/grassroots/TCI) → opex (`opex.py`, raw materials + utilities + Turton
  COM) → profitability (`profitability.py`, NPV/IRR/payback/**LCOP**) → uncertainty
  (`uncertainty.py`, Monte-Carlo P10/P50/P90 + tornado). Validated vs a Turton
  worked example and closed-form finance. `evaluate_economics()` is the cheap
  re-evaluable core (no re-solve) the MC/tornado reuse.
- Examples `05_ammonia_economics` (M3 DoD: costing report + tornado, LCOP
  ~$0.69/kg NH3) and `06_equation_oriented` (M4 DoD: SM vs EO agree to ~1e-10 on
  the recycle flowsheet, then a min-duty-s.t.-recovery optimization).

**M4 note on the EO backend:** built on **scipy** (Newton for the simultaneous
solve, SLSQP for optimization) reusing Caldyr's `unit.solve` + property package,
NOT IDAES/Pyomo/IPOPT (none installed, and a Pyomo model would need a separate
algebraic property model that wouldn't match thermo's PR/NRTL). The scipy backend
keeps both solvers on identical physics — the right call for the cross-check. A
Pyomo/IDAES backend (large-scale, analytic derivatives) is a pluggable future
addition behind the same `backend=` interface.

**Deferred** (additive, not needed for any DoD): Gibbs reactor via Cantera and
kinetic PFR/CSTR reactors (M2); fired-heater-specific costing + feed-effluent heat
integration (M3); Pyomo/IDAES EO backend (M4).

- `api/` (FastAPI bridge, thin transport, no physics): `GET /unit-types` &
  `/property-packages` (palette), `POST /solve` (.flow + backend → resolved
  streams + report), `POST /cost` (full TEA + tornado), `POST /optimize`
  (declarative metric specs), `POST /flow/roundtrip`. Tested headlessly with
  `TestClient`. Run: `python -m uvicorn api.main:app --port 8753`.
- `web/` (Vite + React + TS + `@xyflow/react`, thin client): a flowsheet canvas
  with a unit-op palette, custom nodes (one handle per Port, energy ports amber),
  stream edges, a param panel, Solve → stream table, Cost → economics tab with a
  sensitivity tornado. `canvasToFlow`/`flowToCanvas` keep canvas == `.flow` code.
  Boundary feeds/products are UI-only nodes. Dev port **5273** (API proxied at
  `/api` → 8753). `npm install && npm run build` green; build+solve+cost verified
  in a real browser on the ammonia loop.
  **UI Phase A (2026-06)**: state lives in a Zustand store (`src/store.ts`);
  components split under `src/components/`; undo/redo, copy/paste (Ctrl+C/V/D),
  double-click rename, file save/open + localStorage autosave (UI extras ride in
  the engine-ignored `meta.ui` key), Feed/Product palette + components editor
  (flowsheets buildable from scratch), param validation with units
  (`src/lib/params.ts`), toasts/loading/stale states, Tailwind v4 design tokens
  (`@theme inline`) with light/dark themes. Tests: `npm test` (Vitest round-trip
  for flow.ts) and `npm run test:e2e` (Playwright smoke: solve+cost, undo).

**M5 note:** delivered a working single canvas view (the DoD: build+solve+cost in
the browser). The full BFD/PFD/P&ID three-view toggle and live plots are future
polish, not yet built.

- `caldyr.ai` (the AI layer, **local-first**): an `AgentSession` plus JSON-schema'd
  tools (`tools.py`: new_flowsheet, add_unit, add_feed, connect, solve, cost,
  optimize, stream_table, export_flow, list_*). `dispatch(session, name, args)` runs
  one tool, returns structured data + an LLM-friendly summary, and returns errors as
  `{ok:False,error}` (validated ports, normalized compositions, solve invalidated by
  edits — so models recover gracefully). LLM backends are **pluggable and default to
  local** (`llm.py`): `ollama` (default, no key/cost), `openai` (any
  OpenAI-compatible local server), `anthropic` (opt-in). Select via `run(prompt,
  provider=...)` or `CALDYR_LLM_PROVIDER`/`CALDYR_LLM_MODEL`. `mcp_server.py` exposes
  the tools over MCP for Codex CLI / Claude Desktop
  (`codex mcp add caldyr -- python -m caldyr.ai.mcp_server`). Example `07_ai_agent`
  runs the local agent (or a scripted transcript offline). Verified live: a local
  qwen3 (Ollama) built + solved + costed the ammonia loop end-to-end.

**Do NOT use the Anthropic/Claude API during development — default to local LLMs
(Ollama).** The `anthropic` provider exists but is opt-in only.

**Two real engine fixes landed in M6:** (1) equipment sizing reads duties by the
engine's default id (`unit.port`) when a duty port is unwired, so costing is correct
whether or not a builder connects duty ports (both solvers report unwired duties);
(2) `economics.analyze` raises a clear, actionable error when no product stream is
found (instead of a bare ZeroDivisionError).

**Next — M7: Open-source launch.** Curated component DB, docs site, validated
example library, CI, contribution guide, public repo, first release. The engine,
economics, EO backend, API, web canvas, and AI layer are all in place and green.

## 8. Non-goals (for now)

- Dynamic simulation (steady-state first).
- Rigorous tray-by-tray columns (shortcut Fenske-Underwood-Gilliland first).
- Reproducing Aspen's full component database (start with a small curated set).
- Pretty UI before a correct engine.

## 9. Where things live

- `engine/caldyr/` — the simulation engine (this is the heart).
- `api/` — FastAPI bridge (thin).
- `web/` — React + `@xyflow/react` front end (thin; no physics).
- `examples/` — runnable, validated example flowsheets.
- `tests/` — pytest; each model validated against a cited reference.
- `docs/` — `ARCHITECTURE.md`, `DATA_MODEL.md`, `THERMO.md`, `ECONOMICS.md`.
