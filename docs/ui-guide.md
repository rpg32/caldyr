# Using the app

## Canvas

- **Palette** (left): Feed/Product boundaries and all 36 unit operations
  (served live from the engine's registry — anything registered in
  `caldyr.unitops` appears automatically; see the
  [unit-operations reference](unit-operations.md)). Click to add.
- **Connect** by dragging between ports — inlets on the left of a node, outlets
  on the right; energy ports are amber. Mismatched or occupied ports are refused
  with an explanation.
- **Dynamic ports**: units whose ports derive from params grow handles as you
  edit — a multi-feed `RigorousColumn` (`feeds`) gains `in1..inN`, side draws
  gain `side1..sideN`, a `Balance` grows with `n_inlets`. The canvas refetches
  the port list (`POST /ports`) whenever such params change.
- **Rename** units by double-clicking; streams in the inspector.
- **Views**: BFD (plain blocks, no duty lines), PFD (full detail), P&ID
  (adds an instrumentation overlay synthesized from your Set/Adjust logical ops).
- **Color by** phase (vapor / liquid / two-phase / duty) or temperature
  (blue→red heat map); a legend appears bottom-left.
- **Hover a stream** for a live T/P/flow/phase readout; **double-click** it to
  pin the callout to the canvas.
- **Arrange** auto-lays-out the flowsheet (recycle-aware). **Group** collapses
  the current selection into a block; double-click the block to expand.

## Keyboard

| Keys | Action |
|---|---|
| Ctrl+Z / Ctrl+Y | Undo / redo |
| Ctrl+C / V / D | Copy / paste / duplicate |
| Ctrl+A | Select all |
| Ctrl+S / Ctrl+O | Save / open `.flow` file |
| Shift+drag | Marquee select |
| Del / Backspace | Delete selection |

## Inspector tabs

- **Params** — selected unit's parameters (with units and validation), selected
  stream's state + phase envelope, or — with nothing selected — the flowsheet
  panel: components (with autocomplete), property package, costing product,
  groups, and the **logical ops editor** (Set/Adjust).
  Below a solved unit's params, the **Design results** panel shows the
  engine's `unit.design` output: scalar tables (FUG numbers, fired-heater fuel
  duty, pipe friction breakdown...) and, for staged columns, the converged
  stage profiles as charts — temperature and liquid-composition vs stage.
- **Streams** — stream table (SI / Metric / Field unit sets, CSV export),
  duties, convergence plot, and the mass & energy **balance check**.
- **Econ** — KPIs, installed equipment costs, interactive tornado, Monte-Carlo
  histogram (P10/P50/P90).
- **Opt** — optimization builder over a **physical** objective/constraint
  (duty, stream flow, or component rate) and design variables. Economic
  objectives like LCOP are Python-only (`solver.optimize` with a callable) — see
  [tutorial 3](tutorial-optimization.md).
- **Study** — parameter sweeps with live charts and CSV export.

## Projects

The **Projects** dialog saves named flowsheets in your browser, and offers
templates — from a two-unit starter through the ammonia loop, shortcut
distillation and SMR hydrogen up to three book-validated plants with closed
recycles (Hameed 2025 ch. 15): the **cyclohexane plant** (benzene
hydrogenation, feed-effluent HX, H2 recycle), the **VCM plant** (EDC cracking,
two pressure columns, EDC recycle), and the **DME plant** (methanol
dehydration, two-column train, methanol recycle). The recycle templates ship
with `solver_hints` tear guesses and solve in seconds-to-a-minute. Your
working canvas also autosaves continuously and restores on reload. `.flow`
files saved to disk are plain JSON, engine-compatible, and git-friendly.

## Copilot

The **Copilot** panel chats with a local LLM (Ollama by default) that operates
the same typed tools as the Python API: it can build, edit, solve, cost,
optimize, and explain flowsheets. Every structural edit it proposes is shown as
a **diff against your canvas** — nothing lands until you accept. The
*Explain flowsheet* and *Diagnose solve* buttons work even without an LLM.
