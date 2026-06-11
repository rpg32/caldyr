"""M11 demo: fired heating, air cooling, and heat-integration targeting.

    FEED ─▶ FIRED (fired heater, 300 → 650 K) ─▶ AIRCOOL (650 → 330 K) ─▶ PRODUCT

A methane stream is heated in a direct-fired heater (fuel duty = process duty /
fired efficiency) and then cooled back down against ambient air (no cooling
water — the only utility is fan power). The flowsheet is solved, the equipment
is sized and costed (Turton 4e correlations), and then `pinch_analysis` targets
the heat integration: the heating demand (cold stream, 300→650 K) and the heat
rejection (hot stream, 650→330 K) overlap almost completely, so a feed-effluent
exchanger could displace most of the fired duty — the analysis quantifies
exactly how much, before any exchanger is designed. Run from the repo root:

    python examples/12_heat_integration.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.analysis import pinch_analysis  # noqa: E402
from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.economics import cost_equipment, size_flowsheet  # noqa: E402
from caldyr.thermo import make_package  # noqa: E402
from caldyr.unitops import AirCooler, FiredHeater  # noqa: E402

T_FEED = 300.0      # K
T_HOT = 650.0       # K, fired-heater outlet
T_PRODUCT = 330.0   # K, air-cooler outlet
P = 2.0e6           # 20 bar
FLOW = 100.0        # mol/s methane


def build() -> Flowsheet:
    fs = Flowsheet(components=[Component("methane")], property_package="thermo:PR")
    fs.add(FiredHeater("FIRED", {"T_out": T_HOT, "efficiency": 0.85}))
    fs.add(AirCooler("AIRCOOL", {"T_out": T_PRODUCT}))
    fs.feed("FEED", "FIRED:in1", T=T_FEED, P=P, molar_flow=FLOW, z={"methane": 1.0})
    fs.connect("S1", "FIRED:out", "AIRCOOL:in1")
    fs.connect("PRODUCT", "AIRCOOL:out", None)
    fs.connect("Q_FIRED", "FIRED:duty", None)
    fs.connect("Q_AIRCOOL", "AIRCOOL:duty", None)
    return fs


def main() -> None:
    fs = build()
    report = fs.solve()
    print(f"Converged: {report.converged}")
    q_fired = report.duties["Q_FIRED"]
    q_cool = report.duties["Q_AIRCOOL"]
    fuel = fs.units["FIRED"].design["fuel_duty"]
    print(f"\nFired heater : process duty {q_fired / 1e6:6.2f} MW, "
          f"fuel duty {fuel / 1e6:6.2f} MW (efficiency 0.85)")
    print(f"Air cooler   : heat rejected {-q_cool / 1e6:6.2f} MW "
          f"(fan power {0.02 * abs(q_cool) / 1e3:.0f} kW)")

    # -- equipment sizing + Turton bare-module costing -----------------------
    pp = make_package(fs.property_package, fs.component_ids)
    sizes = size_flowsheet(fs, report, pp)
    print(f"\n{'unit':>8} {'equipment':>14} {'capacity':>16} {'Cbm (2023$)':>12}")
    print("-" * 56)
    for size in sizes:
        cost = cost_equipment(size)
        print(f"{size.unit_id:>8} {size.equipment_type:>14} "
              f"{size.attribute:>10.1f} {size.attribute_name:<5} "
              f"{cost.bare_module:>11,.0f}")

    # -- heat-integration targeting -------------------------------------------
    res = pinch_analysis(fs, report, dt_min=10.0)
    print("\nPinch analysis (dT_min = 10 K):")
    print(f"  current hot utility      : {res.current_hot_utility / 1e6:6.2f} MW")
    print(f"  current cold utility     : {res.current_cold_utility / 1e6:6.2f} MW")
    print(f"  minimum hot utility      : {res.qh_min / 1e6:6.2f} MW")
    print(f"  minimum cold utility     : {res.qc_min / 1e6:6.2f} MW")
    if res.pinch_T_shifted is not None:
        print(f"  pinch temperature        : {res.pinch_T_shifted:.1f} K shifted "
              f"(hot {res.pinch_T_hot:.1f} / cold {res.pinch_T_cold:.1f})")
    else:
        print("  pinch temperature        : none (threshold problem)")
    print(f"  heat-recovery potential  : {res.heat_recovery_potential / 1e6:6.2f} MW "
          f"({100 * res.heat_recovery_potential / res.current_hot_utility:.0f}% of "
          f"the fired duty)")
    print("\n  -> a feed-effluent exchanger recovering this duty would shrink both")
    print("     the fired heater and the air cooler to the targets above.")


if __name__ == "__main__":
    main()
