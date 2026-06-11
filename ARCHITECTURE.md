# Architecture

Caldyr is a layered system. Each layer depends only on the ones below it. The
engine is fully usable headless (Python); the API and web app are thin shells.

```
┌─────────────────────────────────────────────────────────────┐
│  web/  React + @xyflow/react   (canvas, stream tables, plots) │  presentation
├─────────────────────────────────────────────────────────────┤
│  ai/   typed tools over engine  (Claude builds/edits/explains)│  assistance
├─────────────────────────────────────────────────────────────┤
│  api/  FastAPI    (REST/WS; serialize flowsheet ⇄ engine)     │  transport
├─────────────────────────────────────────────────────────────┤
│                        engine/caldyr                         │
│  io        .flow JSON load/save, schema, versioning           │
│  economics sizing → costing → profitability → scaling         │
│  solver    sequential-modular  |  IDAES/Pyomo (equation-oriented)│ orchestration
│  unitops   Mixer, Splitter, Heater, HX, Flash, Pump, Comp,    │
│            Reactor, Column …  (common UnitOp.solve contract)   │
│  thermo    PropertyPackage interface → thermo/CoolProp/Cantera│  physics
│  core      Component, Stream, Port, UnitOp, Flowsheet         │  data model
└─────────────────────────────────────────────────────────────┘
```

## Layer responsibilities

### core — data model
Pure data + connectivity, no physics. `Flowsheet` is a directed graph: nodes are
`UnitOp`s, edges are `Stream`s bound to typed `Port`s. Owns the component list and
the active `PropertyPackage`. See `docs/DATA_MODEL.md`.

### thermo — properties
A single `PropertyPackage` interface so the rest of the engine never imports a
specific backend:

```python
class PropertyPackage(Protocol):
    def enthalpy(self, stream: StreamState) -> float: ...
    def flash_pt(self, stream: StreamState) -> PhaseResult: ...   # P,T flash
    def flash_ph(self, P: float, H: float, z: dict) -> PhaseResult: ...
    def bubble_dew(self, stream: StreamState) -> tuple[float, float]: ...
```

Implementations wrap `thermo`/`chemicals` (general EOS, mixtures), `CoolProp`
(fast pure-fluid/REFPROP-grade), `Cantera` (equilibrium/kinetics). Selecting a
package is per-flowsheet. This is the boundary that keeps us out of the business
of re-deriving thermodynamics.

### unitops — models
Every unit op implements one contract:

```python
class UnitOp(ABC):
    ports: list[Port]
    def solve(self, inlets: dict[str, Stream], pp: PropertyPackage)
        -> dict[str, Stream]: ...        # returns outlet streams
```

This uniform contract is what makes the sequential-modular solver and the AI tool
layer simple. Start with the algebraic ones (Mixer, Splitter, Heater/Cooler,
Valve, Pump, Compressor, Flash drum), then HX, then Reactors (conversion →
equilibrium → Gibbs → kinetic), then shortcut Distillation.

### solver — orchestration
Two backends behind one `solve(flowsheet)` facade:

- **Sequential-modular (M0–M1, our own):** topologically sort the unit ops; for
  cycles, detect tear streams (minimum tear set), then converge recycles with
  direct substitution → Wegstein → Newton fallback. Tractable, transparent, great
  for the AI layer and teaching. This is the default.
- **Equation-oriented (M4+, via IDAES/Pyomo):** assemble the whole flowsheet as
  one nonlinear system and hand it to IPOPT for simultaneous solve and, crucially,
  **optimization and scale-up** (the capability gap Aspen leaves open). Caldyr
  flowsheets translate to IDAES models through an adapter; we do not fork IDAES.

### economics — techno-economics (first-class)
`sizing → costing → profitability → scaling`. Detailed methodology in
`docs/ECONOMICS.md`. Drives both the product (a costing tab in the UI) and the
separate `ideal-system-analysis` project.

### io — persistence
`.flow` = JSON document: schema version, component list, property package,
unit-op nodes (type + params + position), stream edges (endpoints + spec), and
cached solved state. Plain text → git diff/merge/review. Round-trip must be exact.

### api / web / ai
- `api/` (FastAPI): load/solve/save a flowsheet; stream results over WebSocket for
  long solves; no physics.
- `web/` (React + `@xyflow/react`): nodes = unit ops with custom handles per port;
  edges = streams; side panels for params, stream tables, economics, plots.
  Renders three views over one model: **BFD** (blocks), **PFD** (equipment +
  stream table), **P&ID** (instrumentation/control overlay).
- `ai/`: exposes engine actions (`add_unit`, `connect`, `set_param`, `solve`,
  `cost`, `explain_convergence`) as typed tools so Claude operates the engine the
  same way a user does — the "conversational process engineering" differentiator.

## Why this beats the incumbent
Open + free + git-native + scriptable + AI-driven, with TEA and IDAES-grade
optimization in the core. We inherit validated thermo instead of rebuilding it,
so effort concentrates on the experience and the analysis capabilities Aspen
lacks.
