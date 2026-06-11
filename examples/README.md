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
