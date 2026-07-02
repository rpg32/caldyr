# Optimize a flowsheet against an economic objective

**Goal.** Move from *solving* a flowsheet to *optimizing* it: let a solver pick
the design variables that minimize an objective, subject to constraints — first a
**physical** objective (minimum duty at a required recovery), then an
**economic** one (minimum LCOP). Then quantify the answer with **tornado**
sensitivity and **Monte-Carlo** uncertainty. You will do the physical
optimization both in the web **Opt** panel and in Python, and see honestly where
each path stops.

**Time:** ~20 minutes. **You need:** Caldyr installed
([Getting started](getting-started.md)); the web app for the Opt-panel section.
Source: `examples/06_equation_oriented.py` (solve + optimize),
`examples/05_ammonia_economics.py` (tornado + Monte-Carlo),
`examples/13_pyomo_backend.py` (the EO-backend limits). If you have not done the
[ammonia](tutorial-ammonia-loop.md) and
[distillation](tutorial-distillation-tea.md) tutorials, do them first — this one
builds on both.

---

## How optimization works in Caldyr

`caldyr.solver.optimize` wraps scipy's SLSQP around a **full flowsheet solve**:
every time it needs the objective or a constraint, it re-solves the flowsheet at
the current design point and reads a scalar off the result. So an "optimization"
is an outer loop of tens of solves. You give it three things:

- an **objective** — any callable `(fs, report) -> float` to minimize;
- **design variables** — `DesignVar(unit_id, param, lower, upper, initial)`, the
  unit parameters it is allowed to move;
- **constraints** — callables `(fs, report) -> float` that must stay ≥ 0.

Because the objective is an arbitrary Python function, it can be a duty, a
recovery, or a full techno-economic LCOP — anything you can compute from a solved
flowsheet.

---

## Part 1 — a physical objective (min duty at a target recovery)

A flash-with-recycle separates n-pentane from n-octane. We want to choose the
**flash temperature** that **minimizes heating duty** while still recovering at
least 4.2 mol/s of pentane overhead. This is `examples/06_equation_oriented.py`.

### In Python

```python
from caldyr.core import Component, Flowsheet
from caldyr.solver import DesignVar, optimize
from caldyr.unitops import FlashDrum, Mixer, Splitter

def build():
    fs = Flowsheet(components=[Component("n-pentane"), Component("n-octane")],
                   property_package="thermo:PR")
    fs.add(Mixer("MIX", {"dP": 0.0}))
    fs.add(FlashDrum("FL", {"T": 360.0, "P": 101325.0}))
    fs.add(Splitter("SP", {"split": 0.6}))
    fs.feed("FEED", "MIX:in1", T=330.0, P=101325.0, molar_flow=10.0,
            z={"n-pentane": 0.5, "n-octane": 0.5})
    fs.connect("MIXOUT", "MIX:out", "FL:in1")
    fs.connect("VAP", "FL:vapor", None)
    fs.connect("LIQ", "FL:liquid", "SP:in1")
    fs.connect("RECY", "SP:out1", "MIX:in2")
    fs.connect("BOT", "SP:out2", None)
    fs.connect("Q", "FL:duty", None)
    return fs

def c5_overhead(fs):
    v = fs.streams["VAP"]
    return v.molar_flow * v.z["n-pentane"]

fs = build()
res = optimize(
    fs,
    objective=lambda fs, rep: rep.duties["Q"] / 1e3,           # minimize kW
    design_vars=[DesignVar("FL", "T", 340.0, 370.0, initial=360.0)],
    constraints=[lambda fs, rep: c5_overhead(fs) - 4.2],       # ≥ 0  → C5 ≥ 4.2
    solve_kwargs={"tol": 1e-9},
)
print(res.design, res.objective, res.success, res.n_solves)
```

```
baseline T = 360.0 K:  duty = 216.5 kW, C5 overhead = 4.149 mol/s  (infeasible)
optimum   T = 360.8 K:  duty = 222.6 kW, C5 overhead = 4.200 mol/s  (binding)
success=True  flowsheet solves=47
```

The baseline (360 K) is *infeasible* — it only recovers 4.149 mol/s, short of the
4.2 target. The optimizer nudges the flash up to **360.8 K**, which just meets the
recovery (the constraint is **binding** — exactly 4.200), at the minimum duty
that satisfies it (222.6 kW). It took 47 flowsheet solves. The constraint being
active at the optimum is the tell-tale of a correctly posed problem: you are
paying the least duty that still hits spec.

### In the web app (Opt panel)

The **Opt** tab is a builder over exactly this, using declarative *metric specs*
rather than Python callables:

1. Load or build the flash-recycle flowsheet and **Solve** it once.
2. Open the **Opt** tab.
3. **Objective:** sense `min`, metric **duty** of stream `Q`.
4. **Design variable:** unit `FL`, param `T`, bounds `340`–`370`, initial `360`.
5. **Constraint:** metric **component_rate** of `n-pentane` in stream `VAP`,
   operator `>=`, value `4.2`.
6. Press **Optimize**. The panel reports the optimal `FL.T`, the objective, and
   re-solves so the canvas and stream table show the optimum.

<!-- SCREENSHOT: the Opt panel with the duty objective, FL.T design variable, and the n-pentane recovery constraint filled in, showing the optimal result -->

**The honest limit of the panel.** Its objective and constraint metrics are
`duty`, `flow`, and `component_rate` read off the solved flowsheet — the physical
quantities. It cannot take an *economic* objective like LCOP. For that you write
a Python objective, which is Part 2.

---

## Part 2 — an economic objective (min LCOP)

The whole point of building costing into the core is that LCOP is just another
number you can compute from a solved flowsheet — so it can be an *objective*.
Take the ammonia loop from the [first tutorial](tutorial-ammonia-loop.md) and ask:
**what purge fraction minimizes the levelized cost of ammonia?** The purge is a
real economic knob — too much wastes synthesis gas, too little inflates the
recycle and every unit it flows through (the trade-off you swept by hand at the
end of that tutorial).

```python
from caldyr.economics import TEAConfig, analyze
from caldyr.solver import DesignVar, optimize

# build() is the ammonia loop from examples/05_ammonia_economics.py

def lcop(fs, report):
    res = analyze(fs, report, TEAConfig())     # size → cost → LCOP
    return res.profitability.lcop

fs = build()
res = optimize(
    fs,
    objective=lcop,                            # minimize $/kg NH₃
    design_vars=[DesignVar("SPLIT", "split", 0.80, 0.97, initial=0.90)],
    solve_kwargs={"tol": 1e-6, "max_iter": 600},
)
print(res.design, round(res.objective, 4), res.success, res.n_solves)
```

```
optimal split = 0.900   (i.e. a 10% purge)
optimal LCOP  = $0.6895/kg     success=True   n_solves=31
baseline LCOP = $0.6896/kg  (at the template's split = 0.90)
```

The optimizer lands on an **interior optimum at split ≈ 0.900** — a 10% purge,
essentially the template default — and confirms it is a genuine minimum, not a
boundary. That is a useful result: the ammonia template's purge is already
near-optimal for LCOP. But notice the objective barely moves ($0.6895 vs
$0.6896): LCOP is **nearly flat** in the purge fraction, because — as the
[tornado](#part-3-tornado-sensitivity) will confirm — feed price, not loop
configuration, dominates the cost. The optimizer is doing exactly its job and
telling you something true: for this plant, tuning the purge is not where the
money is.

> **Cost matters, so pick your objective.** Each solve of the ammonia loop is
> fast, but `analyze` runs on top of it, so an LCOP optimization is tens of
> solve-plus-cost evaluations (31 here, ~20 s). Keep the design-variable count
> small and the bounds tight; SLSQP is a local method, so a sensible `initial`
> matters.

---

## Part 3 — tornado sensitivity

Optimization finds the best point; **sensitivity** tells you which assumptions
that point actually depends on. `tornado` swings each economic assumption over a
plausible range, one at a time, and measures the effect on LCOP. From
`examples/05_ammonia_economics.py`:

```python
from caldyr.economics import tornado

res = analyze(fs, report, TEAConfig())
for b in tornado(fs, res.sizes, res.config):
    print(f"{b.variable:<16} {b.low_lcop:.3f} -> {b.high_lcop:.3f}  swing {b.swing:.3f}")
```

```
feed price ±15%   $0.618 -> $0.762   (swing $0.144)
capacity factor   $0.699 -> $0.681   (swing $0.017)
capex ±30%        $0.683 -> $0.696   (swing $0.013)
discount rate     $0.688 -> $0.692   (swing $0.004)
```

**Feed price swings LCOP ten times as far as capex does.** This is why the purge
optimization in Part 2 was flat, and it is decision-relevant: for this plant,
negotiating the hydrogen contract matters an order of magnitude more than shaving
equipment cost or tuning the loop. `tornado` and `monte_carlo` re-evaluate the
*cheap* economics core without re-solving the flowsheet, so they are fast even
though optimization is not. In the web app, the **Econ** tab renders this as an
interactive tornado bar chart.

---

## Part 4 — Monte-Carlo uncertainty

A tornado moves one variable at a time. **Monte-Carlo** moves them together —
sampling cost-correlation error and prices jointly — to give a *distribution* of
LCOP rather than a point:

```python
from caldyr.economics import monte_carlo

mc = monte_carlo(fs, res.sizes, res.config, n=3000, seed=1)
print(mc.lcop)     # {'p10': ..., 'p50': ..., 'p90': ...}
```

```
Monte-Carlo LCOP ($/kg, n=3000):  P10 $0.601   P50 $0.691   P90 $0.781
```

Report the band, not the point: the P50 is $0.691/kg (matching the optimized
value), with a P10–P90 spread of $0.60–$0.78/kg. That spread — driven mostly by
the feed-price and correlation-error samples — is the honest uncertainty on the
answer. `seed` makes it reproducible; `n` trades runtime for smoothness. The
**Econ** tab shows this as a histogram with the P10/P50/P90 marked.

<!-- SCREENSHOT: the Econ tab tornado and Monte-Carlo histogram side by side -->

---

## Solver backends, and honest limits on Windows

Optimization and solving both run on a **backend**. Caldyr has three, and the
distinction matters for what works on your machine:

| backend | what it is | availability |
|---|---|---|
| `sequential` (default) | sequential-modular, Wegstein tear | everywhere |
| `equation_oriented` | all residuals at once, scipy nonlinear solve | everywhere |
| `pyomo` | Pyomo/PyNumero grey-box + IPOPT | **needs cyipopt (conda)** |

The first two are pure-Python/scipy and run anywhere Caldyr installs. They agree
to ~1e-9 — `examples/06_equation_oriented.py` solves the same flash-recycle both
ways:

```
backend              iters     VAP n     RECY n    residual
sequential               9   5.54414   6.68379    3.5e-11
equation_oriented       61   5.54414   6.68379    0.0e+00
-> max stream-flow difference: 2.40e-10 mol/s   (no tear stream at all)
```

You can pass `backend="equation_oriented"` to `optimize(...)` too, and the Opt
panel exposes both.

**The Pyomo backend is the honest exception.** `backend="pyomo"` wraps the
flowsheet as a PyNumero `ExternalGreyBoxBlock` solved by IPOPT — but grey-box
models can *only* be solved through the `cyipopt` interface (the plain ASL
`ipopt.exe` from `idaes get-extensions` cannot consume them). And **`pip install
cyipopt` ships no Windows wheels** — it needs the Ipopt headers and an MSVC
build. So on a pip-only Windows box the Pyomo *solve* is unavailable.
`examples/13_pyomo_backend.py` is honest about it and degrades gracefully — it
still constructs and evaluates the model to prove the translation layer works:

```
pyomo backend: live solve unavailable on this machine -
  no grey-box-capable NLP solver is available ... Install cyipopt with
  `conda install -c conda-forge cyipopt` (recommended); pip ships no Windows wheels.
but the Pyomo model constructs and evaluates fine:
  20 variables, 20 grey-box equality constraints (5 streams × [n_C5, n_C8, T, P])
```

To use the Pyomo backend, install cyipopt from conda-forge. **You almost never
need to:** the built-in `sequential` and `equation_oriented` backends cover
solving and optimization on every platform — the Pyomo path exists for
large simultaneous problems and IDAES interop, not for everyday work. Nothing in
these tutorials requires it.

---

## What you learned

- `optimize(fs, objective, design_vars, constraints=...)` runs SLSQP with a full
  solve in the loop; the objective is any callable of the solved flowsheet.
- A **physical** objective (min duty s.t. a recovery constraint) maps cleanly to
  the web **Opt** panel; an **economic** objective (min LCOP) is a Python
  objective that calls `analyze` — the panel's metrics are physical only.
- **Tornado** ranks which assumptions the answer depends on; **Monte-Carlo**
  gives the P10/P50/P90 band. For the ammonia loop both point at feed price — the
  same reason the purge optimum is flat.
- Solve and optimize on `sequential` or `equation_oriented` anywhere; the Pyomo
  backend needs cyipopt from conda and is not required for anything here.

**Back to:** the [tutorials index](tutorials.md) for the full example catalog.
