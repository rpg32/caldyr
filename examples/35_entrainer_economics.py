"""Example 35 — techno-economics (LCOP) of the §9.5.6 anhydrous-ethanol
entrainer plant.

This is Caldyr's wedge: the SAME converged flowsheet that example 34 builds is
handed straight to the techno-economic pipeline
(``economics.analyze(fs, report, TEAConfig)``) — size every unit (the
decanting-condenser column T-100 as a tower + trays + condenser + reboiler, the
water column T-101 likewise), cost them with the Turton bare-module method,
roll up the capital (TCI) and the operating cost (reboiler steam + condenser
cooling + the cyclohexane make-up), and report the **LCOP** (levelized cost of
product, $/kg ethanol) plus a sensitivity tornado and a Monte-Carlo P10/P50/P90
band.

The plant is solved with the proven 30-stage continuation (robust + fast); the
62-stage book-scale column (examples 33/34) sharpens the purity but the economics
— dominated by the reboiler steam and the tower capital — are well represented
here. The decanting condenser is costed through its overhead duty (it cools the
ternary heteroazeotrope from its dew point down to ``condenser_T`` where the two
liquids settle), which the column now exposes to the sizer.

    python examples/35_entrainer_economics.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.economics import (  # noqa: E402
    SizingOptions, TEAConfig, analyze, monte_carlo, tornado,
)
from caldyr.unitops import Makeup, RigorousColumn  # noqa: E402

P_ATM = 101325.0
COMPS = ["ethanol", "water", "cyclohexane"]
# Representative ~2023 merchant prices ($/kg). Ethanol is both the feed (hydrous)
# and the product (anhydrous) — a dehydration plant's value-add is the spec, so
# the LCOP carries the feed-ethanol cost plus the cost of dehydration.
PRICES = {"ethanol": 0.60, "cyclohexane": 1.20, "water": 0.0}


def build() -> Flowsheet:
    fs = Flowsheet(components=[Component(c) for c in COMPS],
                   property_package="thermo:UNIFAC")
    fs.add(RigorousColumn("T100", {
        "n_stages": 30, "feeds": [{"stage": 3}, {"stage": 14}],
        "reflux_ratio": 3.0, "distillate_rate": 4.0,
        "method": "naphtali_sandholm", "reboiled": True,
        "decant_condenser": True, "condenser_T": 305.0,
        "reflux_layer": "organic", "max_iter": 120,
    }))
    fs.add(RigorousColumn("T101", {
        "n_stages": 16, "feed_stage": 8, "reflux_ratio": 1.0,
        "distillate_rate": 2.0, "method": "naphtali_sandholm",
        "reboiled": True, "max_iter": 120,
    }))
    fs.add(Makeup("MK", {"component": "cyclohexane", "target": 8.0,
                         "T": 305.0, "P": P_ATM}))
    fs.feed("FEED", "T100:in2", T=343.0, P=P_ATM, molar_flow=27.78,
            z={"ethanol": 0.87, "water": 0.13})
    fs.connect("ENTR", "MK:out", "T100:in1")
    fs.connect("AQ", "T100:distillate", "T101:in1")
    fs.connect("REC", "T101:distillate", "MK:in1")
    fs.connect("ETOH", "T100:bottoms", None)            # anhydrous-ethanol product
    fs.connect("WATER", "T101:bottoms", None)
    for u in ("T100", "T101"):
        fs.connect(f"{u}_QC", f"{u}:condenser_duty", None)
        fs.connect(f"{u}_QR", f"{u}:reboiler_duty", None)
    fs.solver_hints = {
        "tear_guesses": {
            "ENTR": {"T": 320.0, "P": P_ATM, "molar_flow": 8.0,
                     "z": {"ethanol": 0.25, "water": 0.10, "cyclohexane": 0.65}},
            "REC": {"T": 333.0, "P": P_ATM, "molar_flow": 8.0,
                    "z": {"ethanol": 0.30, "water": 0.10, "cyclohexane": 0.60}},
        },
        "tear_tolerance": 5e-3,
    }
    return fs


def money(x: float) -> str:
    return f"${x/1e6:.2f}M" if abs(x) >= 1e6 else f"${x/1e3:.1f}k"


def main() -> None:
    fs = build()
    t100, t101, mk = fs.units["T100"], fs.units["T101"], fs.units["MK"]
    print("Solving the entrainer plant (continuation)...", flush=True)
    report = None
    for d100, d101, mkt in [(4.0, 2.0, 8.0), (7.0, 4.0, 6.0), (10.0, 6.0, 5.0),
                            (13.0, 9.0, 4.0), (16.0, 12.0, 3.0)]:
        t100.params["distillate_rate"] = d100
        t101.params["distillate_rate"] = d101
        mk.params["target"] = mkt
        report = fs.solve(method="direct", max_iter=30)
    etoh = fs.streams["ETOH"]
    print(f"  T-100 bottoms: x_EtOH={etoh.z['ethanol']:.4f} "
          f"(water {etoh.z['water']:.4f}, cyclohexane {etoh.z['cyclohexane']:.4f})\n",
          flush=True)

    cfg = TEAConfig(
        product_component="ethanol", product_min_fraction=0.85,
        prices_per_kg=PRICES, sizing=SizingOptions(material="CS"),
    )
    res = analyze(fs, report, cfg)

    print("Equipment (Turton bare-module, installed):")
    for size, cost in zip(res.sizes, res.costs):
        attr = f"{size.attribute:.2f} {size.attribute_name.split('_')[0]}"
        print(f"  {size.unit_id:<16} {size.equipment_type:<16} {attr:>14}"
              f" {money(cost.bare_module):>10}")

    cap, op, pr = res.capital, res.opex, res.profitability
    print(f"\nCapital:  ISBL {money(cap.isbl)}  |  grassroots {money(cap.grassroots)}"
          f"  |  TCI {money(cap.tci)}")
    print(f"Opex/yr:  raw materials {money(op.raw_materials)}  |  utilities "
          f"{money(op.utilities)}  |  total {money(op.total)}")
    print(f"Product:  {res.annual_production_kg/1e6:.1f} kt/yr anhydrous ethanol")
    print(f"\n  LCOP = ${pr.lcop:.3f}/kg ethanol   (NPV {money(pr.npv)} @ "
          f"{cfg.discount_rate:.0%}, {cfg.project_years} yr)")

    print("\nLCOP sensitivity (tornado, ±swing $/kg):")
    bars = sorted(tornado(fs, res.sizes, cfg),
                  key=lambda b: abs(b.high_lcop - b.low_lcop), reverse=True)
    for b in bars:
        print(f"  {b.variable:<22} {abs(b.high_lcop - b.low_lcop):.3f}  "
              f"(${b.low_lcop:.3f} – ${b.high_lcop:.3f})")

    mc = monte_carlo(fs, res.sizes, cfg, n=2000, seed=1)
    print(f"\nMonte-Carlo LCOP (n={mc.n}):  P10 ${mc.lcop['p10']:.3f}  "
          f"P50 ${mc.lcop['p50']:.3f}  P90 ${mc.lcop['p90']:.3f}")
    print("\nThe integrated decant column is costed end-to-end (tower + trays + "
          "decanting condenser +\nreboiler); the reboiler steam dominates the "
          "opex — the lever a heat-integration study targets next.")


if __name__ == "__main__":
    main()
