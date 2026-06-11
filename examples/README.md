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
