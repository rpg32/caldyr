# Caldyr

**Open-source, AI-native process simulation & techno-economic analysis.**

Build BFDs / PFDs / P&IDs, run steady-state mass & energy balances, size and
cost equipment, optimize, and explain it all in plain language — from a modern
web UI *and* a scriptable Python API. Think *open Aspen HYSYS + economic
analyzer, but git-native, scriptable, and with an LLM in the loop.*

## Why Caldyr

| | Commercial simulators | Caldyr |
|---|---|---|
| License | Per-seat, expensive | Open source, free |
| Flowsheets | Binary files | Plain JSON (`.flow`) — diff, review, version |
| Automation | GUI-first, COM bridges | API-first: everything is Python |
| Economics | Separate paid product | Sizing → costing → NPV/LCOP/Monte-Carlo in core |
| AI | — | Local-LLM copilot builds, edits, solves & explains |

## What's inside

- **17 unit operations**, each validated against a cited reference: mixers,
  heaters/coolers, splitters, valves, pumps, compressors, expanders, flash
  drums, heat exchangers, three-phase separators, component splitters,
  conversion / equilibrium / Gibbs (Cantera) / CSTR / PFR reactors, and
  shortcut (FUG) distillation.
- **Two solvers, one physics**: sequential-modular (Wegstein-accelerated
  tearing) and equation-oriented (simultaneous), agreeing to ~1e-9.
- **Logical operations**: HYSYS-style Set and Adjust, plus solver hints and
  mass/energy balance diagnostics.
- **Techno-economics**: Turton-correlation costing, CEPCI escalation,
  capital/opex roll-ups, NPV/IRR/LCOP, tornado sensitivity, Monte-Carlo.
- **A polished web canvas**: BFD/PFD/P&ID views, phase/temperature stream
  coloring, auto-layout, groups, optimization & case-study builders, plots,
  and an AI copilot whose every edit arrives as a reviewable diff.

Start with [Getting started](getting-started.md), then the
[app guide](ui-guide.md). Architecture and roadmap live in the repository
(`ARCHITECTURE.md`, `ROADMAP.md`).
