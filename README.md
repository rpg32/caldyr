# Caldyr

**Open-source process simulation with techno-economics built in. Flowsheets that
solve, cost, and optimize — in your browser or in Python.**

Caldyr is a steady-state chemical process simulator — the kind of tool used to
design plants (think Aspen HYSYS) — that is free, git-native, and scriptable. It
solves mass & energy balances over flowsheets, sizes and costs the equipment, and
optimizes, all on validated open thermodynamics.

![A solved ammonia-loop PFD in Caldyr](docs/img/ammonia-loop-solved.png)

> **Want it hosted?** A managed Caldyr — cloud solves, saved projects, collaboration, no setup — is coming.
> [**Join the waitlist →**](WAITLIST_URL) <!-- TODO(launch): replace WAITLIST_URL with the live form link -->

## Why Caldyr

- **Real thermodynamics, validated solvers.** Peng-Robinson / SRK / NRTL / UNIFAC
  via `thermo`, `CoolProp`, and `Cantera` — no invented numbers in the solve path.
  Two backends (sequential-modular + equation-oriented) agree to ~1e-9, and every
  unit op ships with a test against a textbook / DWSIM / Aspen reference.
- **Economics as a first-class citizen.** Equipment sizing → Turton costing →
  capital/opex → NPV/IRR/LCOP → tornado + Monte-Carlo, in the core — the analysis
  Aspen sells as a separate product.
- **Git-native and scriptable.** Flowsheets are plain JSON (`.flow`) that diff and
  review in git; every action in the web canvas is also in the Python API.

## 60-second quickstart

```bash
pip install -e ".[dev,api,ai,kinetics]"      # engine + API + copilot + Cantera
pytest -q                                     # the full validated suite

python -m uvicorn api.main:app --port 8753    # engine API
cd web && npm install && npm run dev          # UI at http://localhost:5273
```

In the UI: **Projects → Templates → Ammonia loop**, press **Solve** (watch live
convergence), then pick a product and press **Cost**. Or headless in Python:

```python
from caldyr.io import load_flow
from caldyr.economics import TEAConfig, analyze

fs = load_flow("flowsheet.flow")
report = fs.solve()                           # or backend="equation_oriented"
tea = analyze(fs, report, TEAConfig(product_component="ammonia"))
print(tea.profitability.lcop)
```

Full tour: `docs/getting-started.md`; architecture in `ARCHITECTURE.md`.

## What's inside

| Area | What you get |
|---|---|
| **Unit operations** | **36 validated models** — mixers, heaters/coolers, fired heaters, air coolers, pumps, compressors, expanders, valves, pipe segments, flash drums & evaporators, heat exchangers (LMTD / ε-NTU) & multi-stream LNG exchangers, three-phase separators & decanters, conversion / equilibrium / Gibbs / CSTR / PFR / Claus reactors, shortcut (FUG) and rigorous (MESH) columns, absorbers/strippers, extraction, gas saturators, solids ops (cyclone, filters), component splitters, and the Balance logical op. |
| **Thermodynamics** | Peng-Robinson & SRK (cubic EOS); NRTL, UNIFAC & UNIFAC-LLE (activity models — azeotropes, VLLE); IAPWS-95 steam — via `thermo` / `CoolProp` / `Cantera`. |
| **Solvers** | Sequential-modular (Wegstein tear convergence) and equation-oriented (simultaneous Newton) on identical physics, agreeing to ~1e-9; SLSQP optimization on top; Set/Adjust logical ops, solver hints, and balance diagnostics. |
| **Economics** | Sizing → Turton 4e bare-module costing (CEPCI escalation) → ISBL/OSBL/TCI capital → raw-materials/utilities/COM opex → NPV/IRR/payback/LCOP → tornado sensitivity + Monte-Carlo P10/P50/P90. |
| **Web app** | Flowsheet canvas with BFD/PFD/P&ID views, phase/temperature stream coloring, auto-layout, groups, undo/copy-paste, a Ctrl+K command palette, projects & templates, optimization & case-study builders, and plots. |
| **AI copilot** *(beta, local-first)* | An optional assistant that builds, edits, and explains flowsheets through typed tools — the engine does all the physics; the model only orchestrates. Runs locally via Ollama by default (no API key); every proposed edit arrives as a reviewable diff. Also exposed over MCP for Claude Desktop / Codex CLI. |

## Validation

Trust is the whole point. **Every one of the 36 worked examples in `examples/`
cites its reference** — textbook, DWSIM, or Aspen — from the ammonia loop and
benzene/toluene columns through crude towers, amine sweetening, Claus sulfur
recovery, and azeotropic entrainer plants. Property methods are checked against
known anchors (NRTL vs the ethanol/water azeotrope; `lnKeq` vs Haber-Bosch
ammonia equilibrium), and the two solver backends cross-check to ~1e-9.
Convergence diagnostics — residual history, tear streams, per-unit mass/energy
closure — are **exposed, not hidden**: see any solve's Streams tab or its
`SolveReport`.

## Contributing & license

Issues and PRs welcome — see [`CONTRIBUTING.md`](CONTRIBUTING.md); every unit op
and property method must ship with a cited validation test. Caldyr is licensed
under the [AGPL-3.0](LICENSE): free for everyone to use, self-host, and modify —
and anyone who offers a modified Caldyr as a network service must share their
changes. Caldyr never vendors third-party GPL code (e.g. DWSIM internals) —
interop happens across its API only.
