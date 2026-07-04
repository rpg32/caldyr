# Caldyr

**Open-source, AI-native process simulation & techno-economic analysis.**

Build BFDs / PFDs / P&IDs, run steady-state mass & energy balances, size and
cost equipment, optimize, and explain it all in plain language — from a modern
web UI *and* a scriptable Python API. Think *open Aspen HYSYS + economic
analyzer, but git-native, scriptable, and with an LLM in the loop.*

!!! tip "Want it hosted?"
    A managed Caldyr — cloud solves, saved projects, collaboration, no setup — is
    coming. [**Join the waitlist →**](WAITLIST_URL)
    <!-- TODO(launch): replace WAITLIST_URL with the live form link -->

## Why Caldyr

| | Commercial simulators | Caldyr |
|---|---|---|
| License | Per-seat, expensive | Open source, free |
| Flowsheets | Binary files | Plain JSON (`.flow`) — diff, review, version |
| Automation | GUI-first, COM bridges | API-first: everything is Python |
| Economics | Separate paid product | Sizing → costing → NPV/LCOP/Monte-Carlo in core |
| AI | — | Local-LLM copilot builds, edits, solves & explains |

## What's inside

- **36 unit operations**, each validated against a cited reference: mixers,
  heaters/coolers, fired heaters, air coolers, splitters, valves, pipe
  segments, pumps, compressors, expanders, flash drums, evaporators, heat
  exchangers, multi-stream (LNG) exchangers, three-phase separators, decanters,
  component splitters, conversion / equilibrium / Gibbs (Cantera) / CSTR / PFR /
  Claus reactors, shortcut (FUG) and rigorous (MESH) distillation, absorbers/
  strippers, reboiled absorbers, liquid-liquid extraction columns, gas
  saturators, solids ops (cyclone, rotary-vacuum filter, baghouse), and the
  HYSYS-style Balance op — see the
  [unit-operations reference](unit-operations.md).
- **Two solvers, one physics**: sequential-modular (Wegstein-accelerated
  tearing) and equation-oriented (simultaneous), agreeing to ~1e-9.
- **Logical operations**: HYSYS-style Set, Adjust and Balance, plus solver
  hints (tear seeding) and mass/energy balance diagnostics.
- **Analysis tools**: pinch targeting, property tables, psychrometrics,
  API 520/526 relief-valve sizing, phase envelopes, optimization, case
  studies — see [Analysis tools](analysis-tools.md).
- **Techno-economics**: Turton-correlation costing, CEPCI escalation,
  capital/opex roll-ups, NPV/IRR/LCOP, tornado sensitivity, Monte-Carlo.
- **A polished web canvas**: BFD/PFD/P&ID views, phase/temperature stream
  coloring, auto-layout, groups, optimization & case-study builders, column
  stage-profile plots, a template gallery of book-validated plants, and an
  AI copilot whose every edit arrives as a reviewable diff.

Start with [Getting started](getting-started.md), then the
[app guide](ui-guide.md). Architecture and roadmap live in the repository
(`ARCHITECTURE.md`, `ROADMAP.md`).
