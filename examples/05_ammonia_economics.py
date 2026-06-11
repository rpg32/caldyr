"""M3 demo: techno-economic analysis of the ammonia synthesis loop.

Takes the solved M2 ammonia loop, sizes and costs every unit (Turton bare-module
method), rolls up capital and operating cost, and reports the levelized cost of
ammonia with a sensitivity tornado and Monte-Carlo P10/P50/P90 bands.

    python examples/05_ammonia_economics.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.economics import TEAConfig, analyze, monte_carlo, tornado  # noqa: E402
from caldyr.unitops import (  # noqa: E402
    EquilibriumReactor, FlashDrum, Heater, Mixer, Splitter,
)

AMMONIA = {"stoich": {"nitrogen": -1, "hydrogen": -3, "ammonia": 2}, "key": "nitrogen"}


def build() -> Flowsheet:
    fs = Flowsheet(
        components=[Component(c) for c in ("nitrogen", "hydrogen", "ammonia", "argon")],
        property_package="thermo:PR")
    fs.add(Mixer("MIX", {"dP": 0.0}))
    fs.add(Heater("PREHEAT", {"T_out": 673.15}))
    fs.add(EquilibriumReactor("RXN", {"reaction": AMMONIA, "T": 673.15}))
    fs.add(Heater("COOL", {"T_out": 250.0}))
    fs.add(FlashDrum("SEP", {"T": 250.0, "P": 2e7}))
    fs.add(Splitter("SPLIT", {"split": 0.90}))
    fs.feed("MAKEUP", "MIX:in1", T=300.0, P=2e7, molar_flow=100.0,
            z={"nitrogen": 0.2475, "hydrogen": 0.7425, "ammonia": 0.0, "argon": 0.01})
    fs.connect("S1", "MIX:out", "PREHEAT:in1")
    fs.connect("S2", "PREHEAT:out", "RXN:in1")
    fs.connect("S3", "RXN:out", "COOL:in1")
    fs.connect("S4", "COOL:out", "SEP:in1")
    fs.connect("PRODUCT", "SEP:liquid", None)
    fs.connect("VAP", "SEP:vapor", "SPLIT:in1")
    fs.connect("RECYCLE", "SPLIT:out1", "MIX:in2")
    fs.connect("PURGE", "SPLIT:out2", None)
    for u in ("PREHEAT", "RXN", "COOL", "SEP"):
        fs.connect(f"Q_{u}", f"{u}:duty", None)
    return fs


def money(x: float) -> str:
    return f"${x:,.0f}"


def tornado_chart(bars, width: int = 30) -> str:
    span = max(b.swing for b in bars) or 1.0
    lines = []
    for b in bars:
        n = max(1, round(b.swing / span * width))
        lines.append(f"  {b.variable:<16} |{'#' * n:<{width}}|  "
                     f"${b.low_lcop:.3f} -> ${b.high_lcop:.3f}  (swing ${b.swing:.3f})")
    return "\n".join(lines)


def main() -> None:
    fs = build()
    report = fs.solve(tol=1e-7, max_iter=400)
    res = analyze(fs, report, TEAConfig())

    print("=" * 64)
    print("  AMMONIA SYNTHESIS LOOP - TECHNO-ECONOMIC ANALYSIS (2023)")
    print("=" * 64)

    print("\nEquipment (Turton bare-module, installed):")
    print(f"  {'unit':<8} {'type':<18} {'size':>12} {'Cbm':>14}")
    for size, cost in zip(res.sizes, res.costs):
        attr = f"{size.attribute:.2f} {size.attribute_name.split('_')[0]}"
        print(f"  {size.unit_id:<8} {size.equipment_type:<18} {attr:>12} {money(cost.bare_module):>14}")

    cap = res.capital
    print("\nCapital:")
    print(f"  ISBL (installed) ........ {money(cap.isbl)}")
    print(f"  OSBL (offsites) ......... {money(cap.osbl)}")
    print(f"  grassroots (fixed cap) .. {money(cap.grassroots)}")
    print(f"  working capital ......... {money(cap.working_capital)}")
    print(f"  total capital (TCI) ..... {money(cap.tci)}")

    op = res.opex
    print(f"\nOperating cost ($/yr, {op.operating_hours:.0f} h/yr):")
    print(f"  raw materials ........... {money(op.raw_materials)}")
    print(f"  utilities ............... {money(op.utilities)}")
    print(f"  fixed (labor/maint/OH) .. {money(op.fixed)}")
    print(f"  total (COM) ............. {money(op.total)}")

    pr = res.profitability
    print("\nProduction & profitability:")
    print(f"  NH3 production .......... {res.annual_production_kg/1e6:.1f} kt/yr")
    print(f"  revenue (@ $0.50/kg) .... {money(res.annual_revenue)}/yr")
    print(f"  LCOP .................... ${pr.lcop:.3f}/kg NH3  (= ${pr.lcop*1000:.0f}/t)")
    irr = f"{pr.irr*100:.1f}%" if pr.irr is not None else "n/a (cost > price)"
    payback = f"{pr.payback_years:.1f} yr" if pr.payback_years else "n/a"
    print(f"  NPV ..................... {money(pr.npv)}   IRR: {irr}   payback: {payback}")

    print("\nLCOP sensitivity (tornado):")
    print(tornado_chart(tornado(fs, res.sizes, res.config)))

    mc = monte_carlo(fs, res.sizes, res.config, n=3000, seed=1)
    print(f"\nMonte-Carlo LCOP ($/kg, n={mc.n}):  "
          f"P10 ${mc.lcop['p10']:.3f}  P50 ${mc.lcop['p50']:.3f}  P90 ${mc.lcop['p90']:.3f}")
    print("\nNote: at merchant H2 = $1.50/kg the toy loop's LCOP exceeds the $0.50/kg")
    print("sale price, so NPV is negative - feed cost dominates, the green-NH3 story.")


if __name__ == "__main__":
    main()
