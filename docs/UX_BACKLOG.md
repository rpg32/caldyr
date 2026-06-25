# Caldyr Web UI — UX backlog (from the 2026-06-25 user "drive" session)

A hands-on walkthrough of the web app (simple heater, then a flash drum) surfaced
the following user-facing gaps. Listed in **recommended fix order**. The engine is
correct throughout — these are almost all **web display / input / transparency**
issues (engine stays SI internally; it already exposes everything needed).

Servers for testing: API `python -m uvicorn api.main:app --port 8753`; web
`cd web && npm run dev` (→ http://localhost:5273, proxies `/api`→8753). The app
autosaves to localStorage and **reloads the last session** — use **New** for a
clean canvas (it pops a native confirm the *user* must dismiss).

---

## 1. (HIGH BUG) Stream composition is not shown ANYWHERE — #54
The #1 output of any separator is invisible in the UI. Confirmed: the Streams
table shows only `T / P / flow / phase`; the selected-stream **side panel** shows
`name/T/P/flow/phase/vapor-frac` (+ a "Phase envelope" button); the double-click
**pinned canvas card** shows `T/P/flow/phase`. None show per-component fractions.
- **Fix:** (a) show full composition (mole frac; ideally mass frac + per-comp flow)
  in the selected-stream side panel; (b) add it to the hover readout / pinned card;
  (c) add per-component columns (expandable rows or horizontal scroll) to the
  Streams table.
- **Where:** `web/` Streams tab component + stream side-panel + canvas stream
  card/hover. Engine has `stream.normalized_z()`; first verify the `/solve` API
  payload already returns per-stream `z` (add if missing).
- **Verify:** flash 50/50 propane/n-butane @ 313.15 K / 700 kPa → vapor
  0.618/0.382, liquid 0.356/0.644 should all be visible.

## 2. (BUG) Numeric input ".5" becomes "05" — #53
Leading-dot decimals are mangled (typing `.5` yields `05`; must type `0.5`).
- **Fix:** numeric field parser should accept `.5`, `0.5`, trailing dot, etc.
- **Where:** `web/` number-input / params input component.

## 3. (FEATURE) Full unit-system support — #49
Inputs are entered in SI (T in K) but Streams default to Field (°F/psia/lbmol/h);
there's a Streams `units: Field/SI` toggle but **inputs ignore it** → inconsistent.
- **Fix:** (a) a user-selectable **default unit system** (SI / Field / custom) that
  applies CONSISTENTLY to every input AND output; (b) **per-field unit override**
  at every value entered — not just dimension swaps (lb/hr↔g/hr) but magnitude
  variants (lb/hr↔lb/s↔lb/min) to match the order of magnitude; (c) show the unit
  next to every value, round-trip correctly (convert only at the I/O edge).
- **Where:** `web/` inputs (`params.ts`/Inspector/Feed form) + result tables
  (Streams/Econ). Engine already SI-internal; the Streams toggle is the seed.
- **STATUS (2026-06-25, merged `f6e7822`):** (a) ✅ + (c) ✅ done — `web/src/lib/
  units.ts` is now a HYSYS-matched dimension registry (T/P/molar+mass flow/power/
  UA); a global unit picker (Toolbar) drives every physical input (params, feed,
  Study/Optimize/Logical bounds via `QuantityInput`/`DimField` + `dimFor`/
  `dimForMetric`) and output (stream table+CSV, side panel, product node, callout,
  duties/balance, pinch, T-profile). Mass flow + mass fractions added (API `/solve`
  returns a per-component `molar_mass` map). **(b) ✅ done (merged `a4c615b`)** —
  per-field unit-override dropdowns on every dimensioned input (params + feed) via
  `QuantityInput`'s optional unit `<select>`; overrides keyed `${nodeId}:${param}`
  ride in `meta.ui` (scoped per-flowsheet, autosaved, travel with the .flow);
  picking the system-default unit clears the override. **#49 fully done.**
  Deferred-as-SI: ΔT_min (interval), relief tool, property-table value outputs.
  Future nicety (not requested): per-column unit override on OUTPUT tables +
  mass-flow ENTRY for feeds.

## 4. (UX) Guided parameter editor for ALL units — #47
Only the column units got a typed param schema; everything else is
**add-parameter-by-name**, so users must already know each unit's param names
(e.g. Heater = `T_out`/`Q`/`dP`, Flash = `T`/`P`). No discoverability.
- **Fix:** extend the typed schema (`api/param_schemas.py` + web `params.ts`) to
  cover every unit op, OR give the "add parameter" box an autocomplete/dropdown of
  that unit's available params with description + unit + default.
- **Source of truth:** engine unit-op `params.get(...)` keys (+ docstrings).
- **STATUS (2026-06-25, merged `210125f`): ✅ DONE.** `api/param_schemas.py` now
  carries a schema for ALL 36 unit ops (was RigorousColumn only), extracted from
  each unit's `params.get` keys/defaults/bounds/docstrings and served by
  `/unit-types`. The web Inspector (`UnitParamsEditor`/`SchemaParamRow`) shows
  every available param with its widget, default, help, unit (+ per-field unit
  picker) and conditional visibility; reset (×) reverts an override; there's no
  free-text add (so typos can't create junk params); unknown/legacy params surface
  under "Extra parameters" with a remove (×); required + out-of-bounds are flagged
  (engine still validates at solve).

## 5. (FEATURE) Expose + make editable ALL cost/heuristic assumptions — #50
Costing uses heuristics and a price catalog the user can't see/edit fully.
- **Fix:** (a) component prices editable per-flowsheet (e.g. nitrogen $0.10/kg from
  `economics/data.py` `PRICES_PER_KG`); (b) transparent + editable OPEX build-up
  (labor, maintenance, overheads, utility prices, capacity factor, discount rate,
  capital factors — the Turton COM); (c) sizing heuristics surfaced/adjustable
  (LMTD, design margins, U values, flooding fraction, HETP); (d) an "assumptions in
  this solve" summary showing exactly which numbers/correlations drove the result,
  with the `data.py` citation.
- **Where:** engine centralizes these in `economics/data.py` + `TEAConfig`
  (already takes `prices_per_kg`); needs API exposure + a web "Assumptions" panel.
- **STATUS (2026-06-25, merged `07b74f0`): ✅ DONE.** Engine: the COM + capital
  factors that were hard-coded in opex.py/capital.py are now `data.CostFactors`
  on `TEAConfig` (+ `utility_prices`); SizingOptions was already on the config;
  defaults reproduce the validated TEA exactly. API: `CostConfig` accepts
  `prices_per_kg` / `utility_prices` / `sizing{}` / `factors{}` overrides and
  `/cost` returns an `assumptions` block (effective config + the prices/utilities
  that applied + sizing + factors + source citations). Web: a collapsible
  **"Assumptions in this solve"** panel in the Econ tab shows + edits every value
  in place (overrides persist per-flowsheet in `meta.ui`), with Re-cost / reset
  and the data.py citations.

## 6. (UX) Jargon glossary + inline term explanations — #52
Lots of acronyms (LCOP, TCI, ISBL/OSBL, CBM, COM, duty, reflux ratio, Murphree,
HETP, LMTD, tear stream, vapor frac…).
- **Fix:** inline tooltips on labels/headers + a searchable glossary/cheatsheet.

## 7. (BUG) Palette hover tooltip truncated — #48
The left-rail unit-op tooltip is cut off (only part of the docstring shows).
- **Fix:** show the full description (wrap/scroll or larger tooltip). `web/` palette.

## 8. (FEATURE) Save configurations + multiple solve cases — #51
Let users save named configs/settings and multiple cases per flowsheet (base vs
cheaper-N2 vs high-capacity), switch/compare (LCOP/NPV/duties/purities side by
side). Persist with the `.flow` (or alongside). Distinct from the single autosave.

---

## After the fixes — re-run these two tutorials to verify

**A. Heater (simple).** Components: `nitrogen`, package thermo:PR.
Feed (T 298.15 K, P 101325 Pa, 10 mol/s, N₂=1) → Heater (`T_out`=473.15 K) →
Product. Solve → S2 = 200 °C, duty ≈ 0.17 MMBtu/h (~50 kW). Cost → LCOP ≈
$0.389/kg, NPV ≈ −$19.8M (cost center; N₂ feed $0.10/kg dominates).

**B. Flash drum.** Components: `propane`, `n-butane`, package thermo:PR.
Feed (313.15 K, 700 kPa, 100 mol/s, 50/50) → Flash (`T`=313.15 K, `P`=700000 Pa)
→ vapor + liquid Products. Solve → vapor 55.06 mol/s (C₃ 0.618), liquid 44.94
mol/s (C₃ 0.356). **This is the case that proves #54** (compositions visible).

**Then:** advanced tutorials — ShortcutColumn (real high-purity C₃/nC₄ split),
reactor + separator + recycle, then costing/optimization on a real plant.

## Notes for the implementer
- CDP **screenshots are flaky on this machine** (Page.captureScreenshot times out);
  `read_page` (accessibility tree) + ref-clicks work. Prefer instruction-driving
  the user, or use `read_page`/console/network tools, not screenshots.
- Web stack: Vite + React + TS + `@xyflow/react`, Zustand store (`web/src/store.ts`),
  components under `web/src/components/`, params/units in `web/src/lib/params.ts`.
  Tests: `npm test` (Vitest), `npm run test:e2e` (Playwright). Keep tsc + vitest +
  `vite build` green; engine ruff + mypy green.
