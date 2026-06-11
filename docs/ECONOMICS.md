# Techno-economic analysis (TEA)

The economics layer turns a *solved* flowsheet into capital cost, operating cost,
and profitability — and supports scaling and uncertainty. This methodology is used
both inside the product (a costing tab) and by the separate `ideal-system-analysis`
project to compare the ideal closed-loop system against the incumbent fossil system.

## Pipeline

```
solved flowsheet
   │  duties, flows, T/P, phases
   ▼
1. SIZING        → equipment dimensions / duty (area, power, volume, stages)
   ▼
2. PURCHASED COST→ Cp from correlations (capacity, material, pressure factors)
   ▼
3. BARE-MODULE   → Cbm = Cp · Fbm   (installation, piping, instrumentation)
   ▼
4. CAPITAL       → ISBL, OSBL, contingency → Total Capital Investment (TCI)
   ▼
5. OPEX          → variable (feed, utilities) + fixed (labor, maintenance)
   ▼
6. PROFITABILITY → NPV, IRR, payback, and levelized cost (LCOP / LCOH / LCOE)
   ▼
7. SCALING + UQ  → six-tenths / learning curves; Monte-Carlo on cost & price
```

## 1. Sizing
From the solved state: heat exchanger area from Q and LMTD/U; compressor/pump
power from the solved work; vessel volume from flow and residence time; column
stages from the shortcut method. Each unit op exposes a `size()` returning a
sizing dict.

## 2–3. Costing correlations
Use published correlations (Turton *Analysis, Synthesis & Design of Chemical
Processes*; Towler & Sinnott). Purchased cost as a function of a capacity
attribute `A`:
```
log10(Cp°) = K1 + K2·log10(A) + K3·(log10(A))²
Cbm = Cp° · (B1 + B2 · Fm · Fp)
```
Escalate to the analysis year with the **CEPCI** index:
`Cost_year = Cost_base · (CEPCI_year / CEPCI_base)`. Keep correlation constants and
CEPCI values in `data/` (cite sources; never hard-code unsourced numbers).

## 4. Capital
- **ISBL** (inside battery limits) = Σ bare-module costs.
- **OSBL** (offsite/utilities) ≈ factor × ISBL.
- Contingency + fees; working capital. → **TCI**.
- For novel/first-of-a-kind plants, add an explicit FOAK factor and track it
  separately from Nth-of-a-kind (NOAK) so the learning-curve story is honest.

## 5. Operating cost
- **Variable:** feedstock, electricity, water, catalysts/solvents — taken directly
  from the solved stream flows × unit prices (prices in `data/`).
- **Fixed:** operating labor, supervision, maintenance (% of ISBL), overhead,
  insurance, property tax.

## 6. Profitability & levelized cost
- Discounted cash flow → **NPV, IRR, payback**.
- **Levelized cost of product** (the headline for the analysis track):
  `LCOP = (annualized TCI + annual OPEX) / annual production`.
  Specializations: LCOH (hydrogen, $/kg), LCOE ($/MWh), $/t-NH₃, $/t-olefin.
- This is the apples-to-apples number to set against incumbent oil/gas routes.

## 7. Scaling & uncertainty
- **Capacity scaling:** `Cost₂ = Cost₁ · (S₂/S₁)^n`, n≈0.6 default (six-tenths),
  per-equipment n where known.
- **Learning curve (deployment scaling):** `Cost(x) = Cost₀ · x^(−b)`,
  `b = −log₂(1−LR)` for learning rate LR (e.g. ~20% for electrolyzers, PV). This
  is the mechanism by which the green-circular system's cost *falls over time* —
  central to the investor story.
- **Monte-Carlo:** sample correlation error, prices, learning rates, capacity
  factor, WACC → distribution of LCOP/NPV, not a single point. Report P10/P50/P90.

## Comparison framework (analysis track)
For each product (H₂, NH₃, olefins, synfuel, steel, freshwater):
1. **Incumbent route** LCOP from fossil feedstock + current capex (well-known).
2. **Caldyr (ideal-system) route** LCOP from the closed-loop flowsheet.
3. Decompose the gap: feedstock vs energy vs capital vs carbon-source cost.
4. Sweep the decisive variables — **electricity price, WACC, learning rate,
   carbon-source cost** — and find the crossover surface (where green ≤ fossil).
5. Express as an investor-facing thesis: capex-heavy, opex-light, cost falling on
   a learning curve vs a fossil incumbent exposed to fuel-price volatility and
   depletion. The destination economics, the transition cost, and the crossover
   conditions — stated explicitly with P10/P50/P90 bands.
