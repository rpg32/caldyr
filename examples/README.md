# Examples
Runnable, validated flowsheets. Each example has a matching test asserting the
balance closes and/or agreeing with a cited reference.

- `01_mixer_heater.py` — M0 architecture spike (acyclic mass+energy balance).
- (M1) flash-with-recycle
- (M2) ammonia synthesis loop  ← feeds the ideal-system-analysis project
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
- `14_absorber.py` — gas absorption: the SO2/water absorber of Hameed (2025)
  sec. 9.1 (sum-rates MESH, stage profiles, Fair-flooding tray design vs the
  book's 1.285 m packed column), a Kremser (book eq. 9.1) cross-check, and
  the book's sec. 9.3.5 reboiled absorber (n-pentane/n-heptane stripping
  tower) reproduced to ~0.1%.
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
