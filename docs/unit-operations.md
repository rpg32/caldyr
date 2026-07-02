# Unit operations

The registered unit operations, grouped by service. Every unit follows the
same contract: typed ports (material or energy), JSON-friendly `params` that
round-trip through `.flow`, a `solve(inlets) -> outlets` implementation on the
flowsheet's property package, and a test validated against a cited reference
(the source column below; full derivations live in each unit's module
docstring under `engine/caldyr/unitops/`).

Conventions that apply everywhere:

- **SI units** in params and streams: K, Pa, mol/s, J, W, m.
- **Duty sign**: energy added to the process stream is positive. Coolers,
  condensers and expanders therefore report negative duties; a Compressor's
  `work` is positive (absorbed), an Expander's negative (delivered).
- **Energy balances close exactly** — stream enthalpies are
  formation-inclusive (see [Thermodynamics](THERMO.md)), so reactor heats and
  adiabatic temperatures come out of plain balances.
- Units that converge an internal calculation (columns, kinetic reactors)
  raise a typed error with diagnostics rather than returning a silent wrong
  answer, and publish their converged design results (stage profiles, FUG
  numbers, fuel duty...) on `unit.design` — shown in the UI's
  **Design results** panel and used by the economics sizers.

## Pressure change & flow

| Unit | Ports | Key params / specs | Validation source | Economics |
|---|---|---|---|---|
| `Mixer` | `in1`, `in2` → `out` | `dP`; outlet PH-flashed at lowest inlet P − dP | First-principles mass/energy closure | No equipment cost |
| `Splitter` | `in1` → `out1`, `out2` | `split` — fraction of flow to `out1`; intensive state unchanged | First principles | No equipment cost |
| `Valve` | `in1` → `out` | `P_out` **or** `dP`; isenthalpic (Joule-Thomson) PH flash | First principles (H conserved) | No equipment cost |
| `PipeSegment` | `in1` → `out` | `length`, `diameter` (required); `roughness`, `elevation_change`, `fittings_K`, `segments`. Single-phase, isothermal; marched pressure profile on `unit.design` | Churchill (1977) all-regime friction; Hameed (2025) ch. 7 pipe network reproduced to ~0.1%; Crane TP-410 / Perry's 8e | Installed pipe cost per L·D^0.74 (Sinnott C&R Vol. 6) |
| `Pump` | `in1` → `out` + `work` | `P_out`, `eta` (default 0.75); incompressible V·ΔP/η model | Turton 4e pump model | Turton, on solved shaft power; electricity opex |
| `Compressor` | `in1` → `out` + `work` | `P_out`, `eta` (default 0.75); adiabatic, isentropic-efficiency path | Isentropic first principles (PS flash) | Turton, on shaft power; electricity opex |
| `Expander` | `in1` → `out` + `work` | `P_out`, `eta` (default 0.80); mirror of Compressor, negative (delivered) work | Isentropic first principles | Turton axial-turbine, on shaft power |

## Heating & cooling

| Unit | Ports | Key params / specs | Validation source | Economics |
|---|---|---|---|---|
| `Heater` | `in1` → `out` + `duty` | `T_out` **or** `Q`, plus `dP`. Heats or cools (Q sign tells) | First principles (Q = n·ΔH) | Exchanger area from Q and the selected hot/cold utility's LMTD; utility opex |
| `FiredHeater` | `in1` → `out` + `duty` | Same spec contract as Heater plus `efficiency` (default 0.85); heats only. `unit.design` carries `fuel_duty = Q/η` | Turton 4e Ch. 8 (fired efficiency 0.80–0.90) | Turton fired-heater on process duty; fired-fuel opex on fuel duty |
| `AirCooler` | `in1` → `out` + `duty` | `T_out` (required), `t_air_in` (default 308.15 K), `approach` (default 10 K, enforced), `dP`, `fan_power_frac` (default 0.02) | GPSA Engineering Data Book; Towler & Sinnott 2e Ch. 12 (approach) | Bare-tube area (U ≈ 40 W/m²K); fan electricity opex — no cooling water |
| `HeatExchanger` | `hot_in`, `cold_in` → `hot_out`, `cold_out` | Exactly one of `duty`, `T_hot_out`, `T_cold_out`, or `UA` (effectiveness-NTU, `arrangement` counterflow/cocurrent); `dP_hot`/`dP_cold`. Outlets from rigorous PH flashes (phase change captured) | Effectiveness-NTU / LMTD textbook forms | Area from Q and LMTD at `hx_overall_U` (default 500 W/m²K); process-process, no utility |
| `MultiStreamExchanger` | `pass{i}_in` → `pass{i}_out` (N≥2 passes, from `passes`) | `passes=[{...}]`, each pass `T_out` **or** `duty` (signed, +=heated) **or** free, plus `dP`; exactly one free pass closes the energy balance. Optional global `min_approach` (MITA, K) or `ua` (W/K) frees a second pass (iterative). `heat_loss` (W). Hot/cold by solution sign, not label | LNG/plate-fin op (Hameed §9.5.2): Weighted (zone) composite-curve analysis with PH-sampled curves — MITA, required UA, LMTD on `unit.design`; a temperature cross raises | Equivalent area `UA/U` (plate-fin costed as shell-and-tube area — order-of-magnitude) |
| `Saturator` | `gas_in`, optional `water_in` → `gas_out`, `liquid_out`, `duty` | `saturant` (default `water`), `relative_humidity` (default 1.0), `T`/`P`/`dP`. Loads the gas with saturant vapour to `RH·y*(T,P)`; an optional pure-saturant `water_in` supply (else auto-supplied); short supply → sub-saturated, surplus → drains | Stream-saturator op (Hameed §10.4): saturation `y*` from a rigorous saturant-rich probe flash; dry gas conserved; isothermal latent-heat duty reported. (Adiabatic-saturation mode = follow-up) | Vertical contacting vessel on gas residence time |

## Separation (single-stage)

| Unit | Ports | Key params / specs | Validation source | Economics |
|---|---|---|---|---|
| `Flash` (FlashDrum) | `in1` → `vapor`, `liquid` + `duty` | `T` + `P` (isothermal PT flash, duty reported) or `P` alone (adiabatic PH flash) | Property-package VLE (see THERMO validation anchors) | Vertical vessel from residence time (default 300 s) |
| `Evaporator` | `in1` → `vapor`, `liquid` + `duty` | `P` plus exactly one of `T`, `duty`, or `vapor_fraction` (bracketed root-find; pure-fluid case handled via saturated enthalpies) | Hameed (2025) §5.2 steam-heated evaporator | Vertical vessel + hot-utility opex on the duty |
| `ThreePhaseSeparator` | `in1` → `vapor`, `liquid_light`, `liquid_heavy` + `duty` | `T` **required** (PT VLLE flash; PH 3-phase deliberately unsupported), `P` (default inlet). PR/SRK only; liquids labeled by mass density | `thermo` FlashVLN; water/n-hexane VLLE structure vs Tsonopoulos (1999) | Horizontal vessel from residence time |
| `ComponentSplitter` | `in1` → `overhead`, `bottoms` + `duty` | `splits` per component, `default_split`, optional `T_/P_overhead/bottoms`. Black-box; the `duty` port keeps energy books honest | None implied (placeholder unit) | Vessel from residence time |

## Columns (multistage)

| Unit | Ports | Key params / specs | Validation source | Economics |
|---|---|---|---|---|
| `ShortcutColumn` | `in1` → `distillate`, `bottoms` + `condenser_duty`, `reboiler_duty` | `light_key`/`heavy_key`, `recovery_light`/`recovery_heavy` (0.99), `rr_factor` (>1, default 1.3), `P`, `partial_condenser`. FUG design results on `unit.design` | Fenske–Underwood–Gilliland–Kirkbride; validated vs Wankat | Tower + sieve trays (Fair flooding diameter) + condenser + reboiler, Turton; utility opex on both duties |
| `RigorousColumn` | `in1` (or `in1..inN` with `feeds`), `distillate`, `bottoms`, optional `side1..sideN` + both duties | `n_stages` (incl. condenser=1, reboiler=N), `feed_stage` or `feeds=[{stage},...]`, `side_draws`, `reflux_ratio`, exactly one of `distillate_rate`/`distillate_to_feed`, `P`, `dP_stage`, `partial_condenser`, `max_iter`. Full stage profiles on `unit.design` | Wang & Henke (1966) bubble-point MESH per Seader 3e ch. 10.3; validated vs the FUG design point | Same as ShortcutColumn (FUG-compatible design keys) |
| `Absorber` | `gas_in` (bottom), `liquid_in` (top) → `vapor_out` (top), `liquid_out` (bottom) | `n_stages` (no condenser/reboiler), `P`, `dP_stage`, `max_iter`. A **stripper is the same unit** with the feeds swapped | Burningham & Otto (1967) sum-rates per Seader 3e ch. 10.4; Kremser closed form (Hameed 2025 eq. 9.1); SO2/water absorber of Hameed sec. 9.1 | Tower + trays (absorption tray efficiency 0.40); no condenser/reboiler |
| `ReboiledAbsorber` | `feed` → `vapor_out`, `bottoms` + `reboiler_duty` | `n_stages` (incl. reboiler), `feed_stage`, exactly one of `vapor_rate`/`boilup_ratio`/`reboiler_duty`, `P`/`P_bottom` or `dP_stage` | Wang-Henke inner loop driven by overhead vapor; Hameed (2025) sec. 9.3.5 stripping tower reproduced to ~0.1% | Tower + trays (stripping efficiency 0.50) + reboiler |
| `ExtractionColumn` | `feed_in` (top), `solvent_in` (bottom) → `extract_out` (top), `raffinate_out` (bottom) + `duty` | `n_stages`, `T` (required, isothermal), `P`, `max_iter`. Stage LLE from `flash_pt_3p` — PR/SRK only; phases tracked by composition continuity | Hameed (2025) sec. 9.4; Seader 3e ch. 8; Kremser cross-check on a PR-quantitative system | Tower + trays from combined liquid throughput (Seader capacity); LL efficiency 0.25 |

## Reactors

| Unit | Ports | Key params / specs | Validation source | Economics |
|---|---|---|---|---|
| `ConversionReactor` | `in1` → `out` + `duty` | `reaction` `{stoich, key}` + `conversion`, or `reactions=[{stoich, key, conversion},...]` applied in order; `T_out` (isothermal) or absent (adiabatic); `dP` | Stoichiometry + formation-inclusive energy balance (closes to machine precision) | Vessel at reactor residence time (default 15 s) |
| `EquilibriumReactor` | `in1` → `out` + `duty` | `reaction` `{stoich, key}`, `T` (isothermal). Single-reaction mass action vs K(T) from `lnKeq` | Haber-Bosch equilibrium vs Smith, Van Ness & Abbott | Vessel |
| `GibbsReactor` | `in1` → `out` + `duty` | `T` (required), `P` or `dP`. Multi-reaction Gibbs minimization (Cantera `gri30.yaml` species mapping); unmappable components with zero flow pass through | Cantera equilibrate; SMR + WGS example flowsheet | Vessel |
| `CSTR` | `in1` → `out` + `duty` | `V` (required), `reactions` (power-law `KineticReaction` dicts), `T` (isothermal) or adiabatic (nested Brent on the energy balance), `dP` | Fogler / Levenspiel closed forms | Vessel on `V` directly |
| `PFR` | `in1` → `out` + `duty` | `V`, `reactions`, optional `T`, `n_steps`, `dP`; stiff LSODA integration; adiabatic T from a PH flash at every volume node | Fogler / Levenspiel closed forms | Vessel on `V` directly |

## Logical

| Unit | Ports | Key params / specs | Validation source | Economics |
|---|---|---|---|---|
| `Balance` | `in1..inN` (param-driven via `n_inlets`) → `out1` (material, or energy in `"heat"` mode) | `mode`: `"mole"`, `"mass"` (needs `z_out`), `"heat"`, `"mole_heat"`, `"mass_heat"`; optional `P` for the recovery flash | Hameed (2025) §6.3 — the HYSYS Balance op, mode for mode | No equipment (logical operation) |

Flowsheet-level **Set** and **Adjust** logical operations (parameter binding
and spec-driven root-finds) are not unit ops — they live in the `.flow`
document under `"logical"` and are described in the
[app guide](ui-guide.md) and `docs/DATA_MODEL.md`.

## Dynamic ports

Most units have a fixed port list, but some derive ports from their params —
a `RigorousColumn` with `feeds=[...]` grows `in1..inN`, with `side_draws`
grows `side1..sideN`; a `Balance` with `n_inlets` grows its inlets. The API's
`POST /ports` returns the per-instance port list for a `{type, params}` pair
(the static `GET /unit-types` palette shows only the defaults), and the web
canvas refreshes a node's handles automatically when its params change.
