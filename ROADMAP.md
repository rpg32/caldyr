# Roadmap

Phased so the engine is always correct and runnable. Each milestone has a
demonstrable, tested deliverable. Don't start a milestone before the previous is
green.

### M0 — Architecture spike (acyclic steady state)
Prove the spine end-to-end on the simplest flowsheet.
- `core`: Component, Stream, Port, UnitOp, Flowsheet.
- `thermo`: one working PropertyPackage (thermo/CoolProp wrapper) — enthalpy + PT flash.
- `unitops`: Mixer, Heater.
- `solver`: sequential-modular for **acyclic** graphs (topo sort).
- `io`: `.flow` load/save round-trip.
- **DoD:** `examples/01_mixer_heater.py` solves; test asserts energy balance closes.

### M1 — Recycles & a real separation
- Tear-stream detection + recycle convergence (direct substitution → Wegstein).
- Unit ops: Splitter, Valve, Pump, Compressor, Flash drum.
- **DoD:** a flash-with-recycle example converges; tolerance + iteration reporting.

### M2 — Reactors & heat integration
- Reactors: conversion, equilibrium, Gibbs (via Cantera), simple kinetic (PFR/CSTR).
- Heat exchanger (LMTD + effectiveness-NTU); heater/cooler with utilities.
- **DoD:** an ammonia-synthesis-loop toy flowsheet (relevant to the analysis track).

### M3 — Economics layer
- Equipment sizing from solved duties/flows.
- Costing: purchased + bare-module (Turton/Towler correlations, CEPCI escalation).
- Capex (ISBL/OSBL/total), opex (variable + fixed), profitability (NPV, IRR, payback,
  LCOP/LCOE), scaling (six-tenths), Monte-Carlo uncertainty. See `docs/ECONOMICS.md`.
- **DoD:** costing report for the M2 flowsheet; sensitivity tornado.

### M4 — Equation-oriented backend (IDAES/Pyomo)
- Adapter: Caldyr flowsheet → IDAES/Pyomo model → IPOPT.
- Simultaneous solve + first optimization (e.g., minimize cost s.t. spec).
- **DoD:** same flowsheet solved both ways agrees within tolerance; one optimization runs.

### M5 — API + web canvas (BFD/PFD/P&ID)
- FastAPI bridge; React + `@xyflow/react` editor; three views over one model.
- Stream tables, param panels, economics tab, plots.
- **DoD:** build + solve + cost a flowsheet entirely in the browser.

### M6 — AI layer
- Typed tools over the engine; Claude builds/edits/explains/debugs flowsheets.
- Natural-language → flowsheet; convergence diagnostics in plain language; auto-reporting.
- **DoD:** "build me an ammonia loop and cost it" produces a solved, costed flowsheet.

### M7 — Open-source launch
- Curated component DB, docs site, validated example library, CI, contribution guide.
- Public GitHub repo, permissive license, first release.

## Parallel: analysis track
The `ideal-system-analysis` project consumes Caldyr as it matures (hand-balances
now → Caldyr-balanced flowsheets at M2 → costed/scaled at M3–M4). The two repos
stay separate; analysis depends on Caldyr, never the reverse.
