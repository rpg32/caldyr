# Tutorials

Three step-by-step tutorials that take you from a blank canvas to a solved,
costed, and optimized flowsheet. They are written for engineers who know
HYSYS/Aspen/DWSIM and want to see how Caldyr does the same work — and the one
thing those tools do not: put techno-economics in the same loop. Do them in
order; each builds on the last.

| # | Tutorial | What you build |
|---|---|---|
| 1 | [Your first flowsheet: ammonia synthesis loop](tutorial-ammonia-loop.md) | A Haber-Bosch loop with reactor, flash, recycle and purge — solved (Wegstein tear), then sized and costed to an LCOP with tornado and Monte-Carlo. Both the web template and a Python script. |
| 2 | [Design and cost a distillation column](tutorial-distillation-tea.md) | A benzene/toluene split: shortcut (FUG) design → rigorous (MESH) verification → Turton costing → a reflux-ratio capital-vs-utility trade-off in dollars. |
| 3 | [Optimize a flowsheet against an economic objective](tutorial-optimization.md) | Minimize duty at a recovery constraint (web Opt panel *and* Python), then minimize LCOP; tornado sensitivity and Monte-Carlo uncertainty; solver-backend limits. |

Every command and number in these tutorials is produced by running the code — the
LCOPs, stream tables, convergence counts and cost roll-ups are real engine
output, not illustrations.

---

## Further examples

The three tutorials are drawn from the `examples/` directory, which is the
authoritative, validated feature tour — 36 runnable scripts, each with a matching
test that asserts the balance closes and/or that results agree with a cited
reference. Run any of them from the repo root:

```bash
python examples/04_ammonia_loop.py
```

The tutorials are built from examples **04, 05** (ammonia loop + economics),
**06** (equation-oriented solve + optimization), **08, 11** (shortcut and
rigorous distillation), and **13** (the Pyomo backend). The full catalog:

| File | Topic | Demonstrates |
|---|---|---|
| `01_mixer_heater.py` | Basics | The simplest flowsheet — two feeds, a mixer and a heater; an acyclic mass + energy balance and a stream table. |
| `02_flash_recycle.py` | Recycle | A flash drum with a liquid recycle — the first torn loop. |
| `03_azeotrope_nrtl.py` | Thermodynamics | Why the property package matters: ethanol/water VLE, ideal vs NRTL, the azeotrope. |
| `04_ammonia_loop.py` | **Recycle / reactors** | **Tutorial 1** — Haber-Bosch synthesis loop: equilibrium reactor, chilled flash, recycle + purge, Wegstein tear, atom balances. |
| `05_ammonia_economics.py` | **Techno-economics** | **Tutorial 1** — sizes and costs the ammonia loop (Turton), rolls up capital/opex, LCOP, tornado, Monte-Carlo. |
| `06_equation_oriented.py` | **Solvers / optimization** | **Tutorial 3** — sequential vs equation-oriented (agree to 1e-9), then optimize flash T for min duty at a recovery constraint. |
| `07_ai_agent.py` | AI copilot | Natural language → a solved, costed flowsheet through the typed AI tools (local LLM). |
| `08_distillation.py` | **Distillation (shortcut)** | **Tutorial 2** — Fenske-Underwood-Gilliland design of a benzene/toluene column + Turton costing. |
| `09_logical_ops.py` | Logical ops | HYSYS-style Set and Adjust, plus the mass/energy balance report. |
| `10_reactors.py` | Reactors | Rigorous reactors: a Gibbs (Cantera) steam-methane reformer feeding a kinetic water-gas-shift PFR. |
| `11_rigorous_column.py` | **Distillation (rigorous)** | **Tutorial 2** — the same split re-rated tray-by-tray (MESH, Wang-Henke), stage profiles, side by side with FUG, costed. |
| `12_heat_integration.py` | Heat integration | Fired heater + air cooler, both costed, then pinch-analysis minimum-utility targets. |
| `13_pyomo_backend.py` | **Solvers (EO)** | **Tutorial 3** — the Pyomo/PyNumero grey-box + IPOPT backend and its honest Windows limits. |
| `14_absorber.py` | Absorption | The SO₂/water absorber (Hameed §9.1), a Kremser cross-check, and a reboiled absorber. |
| `15_utilities.py` | Analysis tools | The "book ch. 2 & 6" utilities batch: property tables, psychrometrics and friends. |
| `16_pipe_relief.py` | Piping & safety | A pipe network (Churchill friction on IAPWS-95 water) and an API 520 relief-valve sizing. |
| `17_extraction.py` | Extraction | Liquid-liquid extraction (acetone/water) and an extractive-distillation column (NRTL VLLE). |
| `18_cyclohexane_plant.py` | Full plant | The cyclohexane plant end to end (Hameed ch. 15.1): benzene hydrogenation, feed-effluent HX, H₂ recycle. |
| `19_crude_assay.py` | Petroleum | Crude assay characterization → pseudo-components (TBP/MW/API curves → 38 hypocomponents). |
| `20_crude_tower.py` | Petroleum | A steam-stripped crude atmospheric tower: side draws, pumparounds, non-reboiled MESH. |
| `21_solids.py` | Solids | Solids operations (Hameed ch. 12): cyclone, rotary vacuum filter, baghouse, a kaolin train. |
| `22_amine_sweetening.py` | Gas treating | A DEA absorber over the reactive acid-gas package. |
| `23_amine_plant.py` | Gas treating | The full amine gas-sweetening plant (Hameed §15.3): absorber + regenerator. |
| `24_amine_plant_flowsheet.py` | Gas treating | The §15.3 amine plant as a native, closed-recycle Flowsheet. |
| `25_amine_plant_reboiled.py` | Gas treating | The §15.3 amine plant with a fully reboiled + refluxed regenerator. |
| `26_claus_sulfur_recovery.py` | Reactors | A Claus sulfur-recovery unit (Hameed §10.2.3). |
| `27_reactive_distillation.py` | Distillation | Reactive distillation (Hameed §9.5.3): reaction and separation on the trays. |
| `28_heteroazeotropic_decant.py` | VLLE | A heterogeneous-azeotrope three-phase decanting stage. |
| `29_packed_column.py` | Distillation | Packed-column sizing (Hameed §9.1 SO₂ absorber): HETP/HTU-NTU, flooding. |
| `30_lng_exchanger.py` | Heat exchange | A multi-stream (LNG) heat exchanger (Hameed §9.5.2). |
| `31_saturator.py` | Unit ops | A gas saturator (Hameed §10.4 stream saturator). |
| `32_fired_heater_design.py` | Unit ops | A fired-heater radiant/convective design split (Hameed §4.3). |
| `33_entrainer_decant_column.py` | Azeotropic | Anhydrous ethanol via an integrated decanting-condenser column (UNIFAC VLLE). |
| `34_entrainer_plant.py` | Full plant | The full closed §9.5.6 anhydrous-ethanol entrainer plant with the cyclohexane recycle. |
| `35_entrainer_economics.py` | Techno-economics | The entrainer plant costed end to end: tower + trays + decanting condenser + reboiler, LCOP, tornado, Monte-Carlo. |
| `36_resid_crude_tower.py` | Petroleum | The full resid-bearing crude atmospheric tower (Naphtali-Sandholm). |

For the engine entry points behind these — solvers, unit ops, analysis tools,
and the TEA pipeline — see [Unit operations](unit-operations.md),
[Analysis tools](analysis-tools.md), and [Economics](ECONOMICS.md).
