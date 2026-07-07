# Examples
Runnable, validated flowsheets. Each example has a matching test asserting the
balance closes and/or agreeing with a cited reference.

- `01_mixer_heater.py` — M0 architecture spike (acyclic mass+energy balance).
- `02_flash_recycle.py` — flash drum with a liquid recycle (M1): a
  pentane/octane feed is mixed, flashed at 360 K / 1 atm under PR, and 60% of
  the flash liquid is recycled to the mixer while the rest leaves as bottoms.
  The recycle stream is torn automatically and converged by Wegstein-
  accelerated direct substitution; the overall mole balance closes and the
  flash duty is reported.
- `03_azeotrope_nrtl.py` — why the property package matters: the ethanol/water
  minimum-boiling azeotrope (~89 mol% ethanol, 78.2 C at 1 atm). Prints the
  bubble-temperature curve under both the PR cubic EOS (monotone toward pure —
  no azeotrope) and the NRTL activity model (reproduces it, ~351.3 K at
  x~0.89), showing why a cubic EOS cannot represent the system.
- `04_ammonia_loop.py` — ammonia synthesis (Haber-Bosch) loop (M2): a makeup
  3:1 H2:N2 stream with 1% argon inert is mixed, preheated, run through an
  EquilibriumReactor at 400 C / 200 bar, cooled and flashed; unreacted vapour
  recycles and a splitter purges 10% to bleed the argon. The recycle is torn
  and converged automatically; reports per-pass conversion, recycle ratio and
  product purity, and closes the N/H/Ar atom balances around the whole loop.
- `05_ammonia_economics.py` — techno-economic analysis (M3) of the M2 ammonia
  loop: sizes and costs every unit by the Turton bare-module method, rolls up
  capital (TCI) and operating cost, and reports the levelized cost of ammonia
  with a sensitivity tornado and Monte-Carlo P10/P50/P90 bands. At merchant H2
  prices the feed cost dominates and NPV is negative (the green-NH3 story).
- `06_equation_oriented.py` — equation-oriented solve + optimization (M4).
  Part 1 solves the flash-with-recycle flowsheet both sequential-modular (tear
  + Wegstein) and equation-oriented (all equations at once, no tear stream) and
  shows they agree. Part 2 uses `solver.optimize` to choose the flash
  temperature that minimizes heating duty subject to a pentane-overhead
  recovery constraint.
- `07_ai_agent.py` — the AI layer (M6): natural language to a solved, costed
  flowsheet through a tool-calling agent. The LLM backend is selectable and
  local-first (Ollama by default; OpenAI/Anthropic via env vars); with no
  backend reachable it runs the same tool calls a model would make as a
  scripted transcript, so the tool layer is demonstrable offline. Also
  registers as an MCP server for Codex and other MCP clients.
- `08_distillation.py` — shortcut distillation (Fenske-Underwood-Gilliland) +
  economics (M7): an equimolar benzene/toluene feed split to 99%/98% key
  recoveries in a ShortcutColumn. The FUG design (N_min, R_min, N, feed stage)
  is read off the unit and the tower, trays, condenser and reboiler are
  Turton-costed, with the levelized cost of benzene reported (dominated by the
  feed raw-material cost).
- `09_logical_ops.py` — flowsheet-level logical operations, Set + Adjust (M8),
  the open analogue of HYSYS SET/ADJUST blocks. Two flash drums in series: two
  Sets slave FL1.T to the upstream heater outlet and FL2.T 15 K above it, and
  an Adjust varies the heater outlet (Brent root-find around the full solve)
  until the heavy product carries 4.0 mol/s. Closes with
  `solver.balance_report` auditing per-unit and overall mass/energy closure.
- `10_reactors.py` — M9 rigorous reactors: steam-methane reforming in a
  GibbsReactor (Cantera Gibbs minimization) followed by a kinetic
  water-gas-shift PFR (power-law kinetics).
- `11_rigorous_column.py` — rigorous tray-by-tray (MESH) distillation: the
  same benzene/toluene split designed with the ShortcutColumn (FUG) and
  re-rated with the RigorousColumn (Wang-Henke bubble-point), side by side
  with converged stage profiles and Turton costing.
- `12_heat_integration.py` — fired heating + air cooling + pinch targeting: a
  FiredHeater (fuel duty = process duty / efficiency) feeds an AirCooler, both
  Turton-costed, then `caldyr.analysis.pinch_analysis` computes minimum-utility
  targets and the heat-recovery potential (the feed-effluent exchanger
  opportunity, quantified before designing any exchanger).
- `13_pyomo_backend.py` — the Pyomo grey-box equation-oriented backend: the
  flash-with-recycle system wrapped as a PyNumero ExternalGreyBoxBlock and
  solved with IPOPT via cyipopt, compared against the sequential-modular
  solution. Where no grey-box NLP solver is installed (cyipopt rarely builds on
  pip-only Windows) it prints the backend's install guidance and still shows the
  Pyomo model constructing and evaluating its residual/Jacobian at the warm
  start.
- `14_absorber.py` — gas absorption: the SO2/water absorber of Hameed (2025)
  sec. 9.1 (sum-rates MESH, stage profiles, Fair-flooding tray design vs the
  book's 1.285 m packed column), a Kremser (book eq. 9.1) cross-check, and
  the book's sec. 9.3.5 reboiled absorber (n-pentane/n-heptane stripping
  tower) reproduced to ~0.1%.
- `15_utilities.py` — the "chapter 2 & 6" utilities batch (M12), Hameed
  (2025). Four mini-demos: a steam point on the `coolprop:Water` IAPWS-95
  package (sec. 2.2), an n-pentane density (T, P) property table under PR
  (sec. 2.1.4), a moist-air humidity point at 30 C / 50% RH (sec. 2.4), and an
  Evaporator + Balance mini-flowsheet (sec. 5.2 / 6.3) that 80%-vaporizes water
  at 70 kPa, recombines the products to prove the balances close, and
  sizes/costs the evaporator.
- `16_pipe_relief.py` — piping & safety, Hameed (2025) ch. 7: the sec. 7.1
  pipe network (700 kg/h of water through valve + 1-in/0.5-in cast-iron pipe
  + 2-hp pump, PipeSegment with Churchill friction on IAPWS-95 water;
  reproduces the book's Fig. 7.13 point-C/point-D pressures to ~0.1%) and the
  sec. 7.3 blocked-outlet steam PSV (API 520 critical vapor sizing: required
  orifice 10.13 vs book 10.17 cm^2, same API 526 "K" selection).
- `17_extraction.py` — liquid-liquid extraction & extractive distillation,
  Hameed (2025) ch. 9: the sec. 9.4 acetone/water extractor washed by
  3-methylhexane (ExtractionColumn LLE cascade; reproduces the book's
  bottoms x(water) = 1.0 and prints stage profiles + Seader-throughput
  tower sizing), then the sec. 9.5.5 extractive-distillation structure on a
  two-feed NRTL RigorousColumn — the same column that pins below the
  acetone/methanol azeotrope (0.789 under ChemSep NRTL) reaches 88 mol%
  acetone with water fed above the main feed (the book's phenol solvent has
  no ChemSep NRTL parameters; substitution documented in the example).
- `18_cyclohexane_plant.py` — the full cyclohexane production plant, Hameed
  (2025) ch. 15.1: benzene hydrogenation (98% conversion, adiabatic) with a
  feed-effluent heat exchanger, cooler, letdown flash, 70% vapour recycle with
  recompression, a degassing flash, and an 11-stage RigorousColumn recovering
  cyclohexane at the bottom (book: x ~ 0.998 at ~184 C / 1450 kPa). Two
  interlocking recycles — the H2 loop and the FEHE thermal loop — stress the
  tear solver on a real plant topology; the torn streams are seeded with
  engineering estimates via `solver_hints`.
- `19_crude_assay.py` — petroleum assay characterization (M14), Hameed (2025)
  ch. 10 "Refinery Process": the sec. 10.2.1 worked crude assay (TBP + MW +
  API curves, i-C4/n-C4/i-C5/n-C5 light ends) characterized by
  `caldyr.assay.characterize_assay` into the book's 38 hypocomponents
  (28 x 25 F + 8 x 50 F + 2 x 100 F cut scheme), each with NBP/SG/MW and
  Kesler-Lee Tc/Pc/omega/Cp constants (citations in `caldyr/assay.py`); the
  pseudo-components then run as ordinary PR species through the book's column
  front end (450 F crude -> 650 F furnace -> 65 psia flash) with a monotone
  light/heavy split and the resid staying in the liquid.
- `20_crude_tower.py` — steam-stripped crude atmospheric tower (M15), after
  Hameed (2025) sec. 10.2.1: preflash drum -> furnace -> a partial-condenser,
  NO-reboiler RigorousColumn with open stripping steam, kerosene/diesel liquid
  side draws and pumparounds (as `stage_duties`). Runs the non-reboiled
  bubble-point MESH (formation-offset-conditioned envelope energy balances, the
  open-steam bottom damped out of its limit cycle); the reflux ratio is an
  OUTPUT, the mass balance closes machine-exact, and the naphtha cut lands
  within ~3% of the assay's TBP-implied yield. Honestly scoped: a *light*
  synthetic crude (the cubic-EOS `thermo` backend returns unphysical hV < hL
  bubble states for resid-range pseudo-components, which no MESH method can
  resolve), side draws not side strippers, and `method='sum_rates'` limit-
  cycles (documented in the example/test).
- `21_solids.py` — solids operations (M15), Hameed (2025) ch. 12 "Solid
  Operation": the sec. 12.1 cyclone (Lapple cut-diameter + grade efficiency at
  the book's Stairmand "High Output" geometry: 98.4% vs the book's 95% design
  spec, with Shepherd-Lapple dP), the sec. 12.2 rotary vacuum filter (McCabe
  continuous cake filtration; reproduces the book's 2.4368 m^2 area and
  0.96957 m drum width with the backed-out cake resistance), the sec. 12.3
  baghouse (air-to-cloth sizing + filter-drag filtration time to a 2 kPa
  dirty-bag limit), and the book's p. 409 cyclone+baghouse kaolin train
  (500 kg/h past the 75% cyclone, ~100% final capture). Solids ride as
  ordinary components; the PSD is a *unit* param (v1 particle model — see
  `caldyr/unitops/solids.py`).
- `22_amine_sweetening.py` — amine gas sweetening (M16), Hameed (2025)
  sec. 15.3: a DEA absorber scrubbing CO2 and H2S from a sour natural gas over
  the modified Kent-Eisenberg reactive acid-gas package (`amine:DEA`, the open
  analogue of HYSYS's "Acid Gas - Chemical Solvents"). A plain equilibrium-
  stage absorber solved by sum-rates (Burningham-Otto); reports acid-gas
  removal, rich loading and the heat-of-absorption temperature bulge, then a
  second MDEA case showing H2S-selectivity appears only once a Murphree
  efficiency limits CO2 mass transfer.
- `23_amine_plant.py` — the full amine gas-sweetening plant (M16), sec. 15.3:
  absorber + steam-stripping regenerator + lean-amine recycle. The regenerator
  is solved by Naphtali-Sandholm simultaneous correction (sum-rates limit-
  cycles and the bubble-point method degenerates on reactive desorption); the
  lean recycle is converged by damped direct substitution with a water makeup.
  Part 2 adds a partial reflux condenser to the regenerator top stage to dry
  the acid-gas product, contrasting the overhead water carried off with and
  without it.
- `24_amine_plant_flowsheet.py` — the sec. 15.3 amine plant as a native
  Flowsheet (M16): the same loop as example 23, but every block is a unit op
  and `fs.solve()` tears and converges the lean-amine recycle automatically. An
  in-loop `Makeup` controller tops the circulating water back to target
  analytically each sweep (the make-up rate an outer Adjust would search for),
  which pins the otherwise-unstable water inventory of an open-steam stripper so
  plain direct substitution converges; the overall component balance is
  verified.
- `25_amine_plant_reboiled.py` — the sec. 15.3 amine plant with a fully
  reboiled + refluxed regenerator (M16): a reflux condenser (dry acid-gas
  product) plus a reboiler duty replacing the open stripping steam, so the
  water loop closes and the make-up collapses to ~1 mol/s (vs ~84 with open
  steam). Pinning both end temperatures is FD-Jacobian-stiff, so the
  regenerator uses `reboiler_duty` with an internal open-steam-proxy warm start;
  the lean recycle is torn by the sequential solver with the `Makeup`
  controller as in example 24.
- `26_claus_sulfur_recovery.py` — a Claus sulfur-recovery unit (M17), Hameed
  (2025) sec. 10.2.3: the H2S-rich acid gas from amine sweetening (examples
  23-25) fed to an adiabatic thermal furnace and two catalytic converters, each
  followed by a SulfurCondenser. The sulfur allotropes (S2 above ~800 K, S8
  below) ride the NASA ideal-gas package (`nasa:claus`) since the cubic EOS
  lacks S2 critical constants; reactor equilibria come from Cantera. Setting the
  air to oxidise exactly a third of the H2S drives the 2:1 H2S:SO2 optimum and
  ~98% sulfur recovery, closing the sulfur-atom and energy balances.
- `27_reactive_distillation.py` — reactive distillation (M18), Hameed (2025)
  sec. 9.5.3: `RigorousColumn` running reversible kinetic reactions on a band
  of trays (per-stage generation enters the MESH balances, heat of reaction
  carried by the formation-inclusive stage enthalpies). Part 1 is the book's
  methyl-acetate synthesis (methanol + acetic acid <=> methyl acetate + water
  on `thermo:NRTL`, the book's Wilson analogue); Part 2 is toluene
  disproportionation on `thermo:PR` as an all-aromatic cross-check reaching
  ~60% conversion. Each is run with and without tray holdup to isolate the
  conversion the reaction adds.
- `28_heteroazeotropic_decant.py` — heterogeneous-azeotrope (three-phase)
  decanting: a condensed water/n-butanol heteroazeotrope overhead settled into
  an entrainer-rich organic layer (refluxed) and an aqueous product layer by the
  `Decanter` unit's VLLE flash. Demonstrated on the PR cubic EOS (robust
  water/organic structure) and on a new NRTL isoactivity liquid-liquid flash
  with illustrative LLE parameters (organic ~0.46 / aqueous ~0.99 water, near
  the experimental 0.51/0.98), giving the activity package a three-phase
  capability.
- `29_packed_column.py` — packed-column sizing, Hameed (2025) sec. 9.1 SO2
  absorber: setting `internals='packed'` on any column sizes it from the
  Eckert/Strigle generalized pressure-drop correlation (flood-limited diameter)
  and HETP rules for bed height, and the economics layer costs it as a tower
  shell + random-packing bed (Perry 8e Tables 14-13 / 14-18) instead of trays.
  Sizes the book's SO2/air/water absorber both ways; the packed (50 mm metal
  Pall at 70% flood) and tray (Fair at 80% flood) diameters bracket the book's
  1.285 m to ~15%.
- `30_lng_exchanger.py` — multi-stream (LNG) heat exchanger, Hameed (2025)
  sec. 9.5.2: a single plate-fin `MultiStreamExchanger` core carrying several
  passes, each with an outlet-temperature or duty spec; whether a pass runs hot
  or cold falls out of the solution. The Weighted (zone) method builds the
  hot/cold composite curves with rigorous PH flashes (capturing phase change)
  and reports the minimum internal approach (MITA) and required UA. A second run
  adds a minimum-approach spec, freeing a second pass and solving iteratively to
  hold the curves 5 K apart.
- `31_saturator.py` — gas saturator, Hameed (2025) sec. 10.4 Stream Saturator:
  the `Saturator` loads a gas with a condensable up to its saturation partial
  pressure at the operating T/P (or a chosen relative humidity), with the latent
  heat reported on the duty port. Saturates dry nitrogen with n-hexane at 300 K
  / 1 atm (PR predicts the hexane vapour pressure to within a few percent —
  water VLE under PR is the documented caveat), then shows the relative-humidity
  knob and an excess-liquid case where the surplus drains.
- `32_fired_heater_design.py` — fired-heater radiant/convective design split,
  Hameed (2025) sec. 4.3: the sec. 4.3.4 problem heating a 1250 kmol/h C1-C5
  feed from 50 to 300 C in a methane-fired heater at 100% excess air and 60%
  fired efficiency. Adding a fuel/excess-air spec to `FiredHeater` runs the
  firebox combustion and the radiant/convective split (Hameed eqs. 4.8-4.9;
  Lobo-Evans furnace balance): fuel and air rates (~62 kmol/h fuel, matching the
  book), flue-gas composition, flame/bridgewall/stack temperatures, the
  radiant-vs-convective duty division and the tube areas.
- `33_entrainer_decant_column.py` — anhydrous ethanol via an integrated
  decanting-condenser column at book scale (62 stages), Hameed (2025)
  sec. 9.5.6: a cyclohexane-entrainer heteroazeotropic column whose ternary
  overhead is settled inside the unit (`RigorousColumn` `decant_condenser=True`),
  the organic layer refluxed in full and the aqueous layer drawn — keeping the
  entrainer circulating so the Naphtali-Sandholm solve stays tractable. Uses the
  predictive UNIFAC VLLE package and two continuations (a cheap 30-stage cold
  solve `warm_start_from`'d onto 62 stages, then a distillate-rate ramp) to
  reach ~99% ethanol bottoms.
- `34_entrainer_plant.py` — the full closed sec. 9.5.6 anhydrous-ethanol
  entrainer plant: the integrated decanting-condenser column T-100 plus a
  water-recovery column T-101 and two recycles (the internal organic reflux and
  the ethanol+cyclohexane distillate returned to T-100). Grown to book scale by
  flowsheet continuation, stage-count continuation (30 -> 62 stages in place)
  and cyclohexane make-up inventory control. Honestly scoped: the closed loop
  reaches ~0.91 ethanol / ~0.53 water here; full book parity sits past a
  Naphtali-Sandholm high-draw turning-point fold that the open-loop column
  (example 33) clears but the recycle does not in one solve. SLOW — not part of
  the fast test suite.
- `35_entrainer_economics.py` — techno-economics (LCOP) of the sec. 9.5.6
  entrainer plant: the same converged flowsheet as example 34 handed to
  `economics.analyze` — sizing both columns (tower + trays + decanting condenser
  + reboiler), Turton-costing them, rolling up capital and operating cost
  (reboiler steam + condenser cooling + cyclohexane make-up), and reporting the
  levelized cost of anhydrous ethanol with a sensitivity tornado and a
  Monte-Carlo P10/P50/P90 band. Solved on the fast 30-stage continuation; the
  reboiler steam dominates the opex.
- `36_resid_crude_tower.py` — the flagship resid-bearing crude atmospheric
  tower (M15), Hameed (2025) sec. 10.2: a steam-stripped atmospheric column on a
  heavy crude whose TBP curve reaches a ~1050 F resid endpoint — preflash drum,
  furnace on the flash liquid, recombined feed into a partial-condenser column
  with no reboiler, open stripping steam, kerosene/diesel liquid side draws and
  pumparounds (`stage_duties`). The non-volatile resid has no bubble point at
  tower pressure so the bubble-point method fails structurally; Naphtali-Sandholm
  (all MESH equations at once, warm-started by inside-out) closes it to machine
  precision, with the naphtha cut near the assay's TBP-implied yield.
</content>
</invoke>
