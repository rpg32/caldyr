# Analysis tools

Beyond solving and costing, Caldyr ships a set of analysis tools — some
operate on a solved flowsheet, some are standalone engineering calculations.
For each: what it does, the engine entry point, the API endpoint (if exposed),
and where it appears in the web app.

| Tool | Engine | API | Web UI |
|---|---|---|---|
| Pinch analysis | `caldyr.analysis.pinch_analysis` | — | — (Python only) |
| Property tables | `caldyr.analysis.property_table` | — | — (Python only) |
| Psychrometrics | `caldyr.analysis.humidity` | — | — (Python only) |
| Relief-valve sizing | `caldyr.analysis.relief_vapor` / `relief_liquid` | — | — (Python only) |
| Balance report | `caldyr.solver.balance_report` | `POST /balance` | Streams tab → balance check |
| Phase envelopes | `PropertyPackage.bubble_dew` | `POST /envelope` | Stream inspector → *Phase envelope* |
| Optimization | `caldyr.solver.optimize` | `POST /optimize` | Opt tab |
| Case studies | (scripted solves) | `POST /solve` per point | Study tab |
| Monte-Carlo & tornado | `caldyr.economics.monte_carlo` / `tornado` | `POST /cost` | Econ tab |
| AI diagnostics | `caldyr.ai` tools | `POST /ai/tool` | Copilot buttons |

## Pinch analysis (heat-integration targeting)

The problem-table algorithm of Linnhoff & Flower (per Kemp 2e; Smith 2005):
given hot/cold streams and a minimum approach `dt_min`, compute minimum hot-
and cold-utility targets, the pinch temperature, composite curves, and the
heat-recovery potential — *before* designing any exchanger network.

- `pinch_analysis(fs, report, dt_min=10.0)` extracts the thermal streams from
  a **solved** flowsheet (every duty-carrying unit contributes one stream;
  column condensers/reboilers contribute from `unit.design`; near-isothermal
  duties are widened to 1 K segments; process-process HeatExchangers count as
  already integrated, so the targets are over the *remaining* problem).
- `pinch_from_streams(streams, dt_min=10.0)` targets a plain list of
  `{"Tin", "Tout", "Q"}` specs with no flowsheet at all.

Validated against Kemp 2e. See `examples/12_heat_integration.py`.

## Property tables

`property_table(pp, z, T=..., P=..., props=...)` — stream properties over a
(T, P) grid, the HYSYS *Stream Analysis → Property Table* equivalent
(Hameed 2025 §2.1.4): plot-ready 2-D arrays of density, enthalpy, vapor
fraction, etc. for any property package and composition. Failed flash points
become NaN (logged) instead of killing the grid. See `examples/15_utilities.py`.

## Psychrometrics (humidity)

`humidity(T, P, rh=... | w=... | t_wb=... | t_dp=...)` — the moist-air state
from any one humidity specification (relative humidity, humidity ratio,
wet-bulb, or dew point), per Hameed 2025 §2.4. Backed by CoolProp's
`HAPropsSI` (ASHRAE RP-1485 real-gas humid-air model — psychrometric-chart
quality), independent of the flowsheet property package.

## Relief-valve sizing

API Standard 520 Part I / API 526 PSV orifice sizing, the HYSYS *Safety
Analysis → PSV sizing* equivalent (Hameed 2025 sec. 7.3):

- `relief_vapor(W, T, M, Z, k, P1, ...)` — critical (choked) vapor flow,
  identical to API 520's SI equation to ~0.02%; subcritical backpressure
  raises a typed error rather than mis-sizing.
- `relief_liquid(W, rho, P1, P2, ...)` — capacity-certified liquid valves.

Both return the required effective discharge area and the standard API 526
orifice letter (D through T). Validated against the book's worked
blocked-outlet steam PSV (10.13 vs 10.17 cm², same orifice selection); see
`examples/16_pipe_relief.py`.

## Balance report

`caldyr.solver.balance_report(fs)` — per-unit and overall mass & energy
closure of a solved flowsheet, the first thing to look at when results smell
wrong. Exposed at `POST /balance` (solves, then reports) and shown in the web
app's **Streams** tab as the balance check.

## Phase envelopes

Bubble/dew temperature vs pressure for any stream's composition, computed by
looping the property package's `bubble_dew` over a log-spaced pressure grid.
`POST /envelope` (`{flow, stream, n, p_min?, p_max?}`); in the UI, select a
stream and press **Phase envelope** in the inspector. Points where the flash
fails (e.g. supercritical) are skipped.

## Optimization

`caldyr.solver.optimize(fs, objective, design_vars, constraints=...)` —
minimize any callable of the solved flowsheet over unit parameters (scipy
SLSQP with a full solve in the loop, auto-scaled objective, bounds per
`DesignVar`). `POST /optimize` takes declarative metric specs (`duty`,
`flow`, `component_rate` of a stream) for the objective and constraints; the
**Opt** tab is a builder over the same specs. See
`examples/06_equation_oriented.py`.

## Case studies (parameter sweeps)

The **Study** tab sweeps any unit parameter over a range, re-solving each
point through `POST /solve`, with live charts and CSV export. Headless, the
same thing is a Python loop over `fs.units[...].params[...] = x; fs.solve()`.

## Monte-Carlo & tornado (economics)

Part of the TEA pipeline (see [Economics](ECONOMICS.md)):
`caldyr.economics.tornado(fs, sizes, cfg)` produces one-at-a-time sensitivity
bars on LCOP; `monte_carlo(fs, sizes, cfg, n=2000)` samples cost-correlation
error and prices for P10/P50/P90 LCOP/NPV distributions. Both re-evaluate the
cheap economics core without re-solving the flowsheet. Exposed on
`POST /cost` (`tornado: true`, `monte_carlo: n`); the **Econ** tab shows the
interactive tornado and the Monte-Carlo histogram.

## AI diagnostics

`describe_flowsheet` and `explain_convergence` (in `caldyr.ai`) — plain-text
narrations of the flowsheet structure and of the last solve's convergence
behaviour. They run without any LLM: `POST /ai/tool` dispatches one tool
statelessly, and the Copilot panel's *Explain flowsheet* / *Diagnose solve*
buttons call exactly that. See [AI copilot & MCP](ai.md).
