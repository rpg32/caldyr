"""M7 demo: shortcut distillation (Fenske-Underwood-Gilliland) + economics.

A benzene/toluene splitter: an equimolar, near-saturated-liquid feed at 1 atm
is separated to 99%/98% key recoveries in a ShortcutColumn. The FUG design
(N_min, R_min, N, feed stage) comes straight off the unit's `design` attribute;
the tower, trays, condenser and reboiler are then sized and costed (Turton
bare-module method) and the levelized cost of benzene reported.

    python examples/08_distillation.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.economics import TEAConfig, analyze  # noqa: E402
from caldyr.unitops import ShortcutColumn  # noqa: E402

P_ATM = 101325.0


def build() -> Flowsheet:
    fs = Flowsheet(components=[Component("benzene"), Component("toluene")],
                   property_package="thermo:PR")
    fs.add(ShortcutColumn("COL", {
        "light_key": "benzene", "heavy_key": "toluene",
        "recovery_light": 0.99, "recovery_heavy": 0.98,
        "rr_factor": 1.3, "P": P_ATM,
    }))
    fs.feed("FEED", "COL:in1", T=365.0, P=P_ATM, molar_flow=100.0,
            z={"benzene": 0.5, "toluene": 0.5})
    fs.connect("DIST", "COL:distillate", None)
    fs.connect("BOT", "COL:bottoms", None)
    fs.connect("QC", "COL:condenser_duty", None)
    fs.connect("QR", "COL:reboiler_duty", None)
    return fs


def money(x: float) -> str:
    return f"${x:,.0f}"


def main() -> None:
    fs = build()
    report = fs.solve()
    assert report.converged
    col = fs.units["COL"]
    d = col.design

    print("=" * 64)
    print("  BENZENE/TOLUENE SHORTCUT COLUMN (Fenske-Underwood-Gilliland)")
    print("=" * 64)

    print("\nDesign (FUG):")
    print(f"  relative volatility ..... alpha_LK = {d['alpha']['benzene']:.3f} "
          f"(PR, geometric mean top/bottom)")
    print(f"  feed quality ............ q = {d['q']:.3f}")
    print(f"  minimum stages .......... N_min = {d['N_min']:.2f}   (Fenske)")
    print(f"  minimum reflux .......... R_min = {d['R_min']:.3f}  (Underwood, "
          f"theta = {d['theta']:.4f})")
    print(f"  operating reflux ........ R = {d['R']:.3f}  (= 1.3 R_min)")
    print(f"  theoretical stages ...... N = {d['N']:.1f}    (Gilliland-Molokanov)")
    print(f"  feed stage from top ..... {d['feed_stage']}      (Kirkbride)")

    dist, bot = fs.streams["DIST"], fs.streams["BOT"]
    print("\nProducts (at column pressure, 1 atm):")
    print(f"  {'stream':<12} {'mol/s':>8} {'T (K)':>8} {'x_benzene':>10} {'x_toluene':>10}")
    for name, s in (("distillate", dist), ("bottoms", bot)):
        print(f"  {name:<12} {s.molar_flow:>8.1f} {s.T:>8.1f} "
              f"{s.z['benzene']:>10.4f} {s.z['toluene']:>10.4f}")

    print("\nDuties:")
    print(f"  condenser ............... {report.duties['QC'] / 1e6:8.2f} MW  (heat removed)")
    print(f"  reboiler ................ {report.duties['QR'] / 1e6:8.2f} MW  (heat added)")

    res = analyze(fs, report, TEAConfig(product_component="benzene",
                                        product_min_fraction=0.9))
    print("\nEquipment (Turton bare-module, installed, 2023):")
    print(f"  {'item':<16} {'type':<18} {'size':>14} {'Cbm':>12}")
    for size, cost in zip(res.sizes, res.costs):
        attr = f"{size.attribute:.1f} {size.attribute_name.split('_')[0]}"
        if size.quantity > 1:
            attr = f"{size.quantity} x {attr}"
        print(f"  {size.unit_id:<16} {size.equipment_type:<18} {attr:>14} "
              f"{money(cost.bare_module):>12}")

    op = res.opex
    print(f"\nCapital (TCI) ............. {money(res.capital.tci)}")
    print(f"Opex (COM, $/yr) .......... {money(op.total)}  "
          f"(utilities {money(op.utilities)})")
    print(f"Benzene production ........ {res.annual_production_kg / 1e6:.1f} kt/yr")
    print(f"LCOP ...................... ${res.profitability.lcop:.3f}/kg benzene")
    print("\nNote: LCOP here is dominated by the raw-material cost of the feed -")
    print("the column itself (capital + utilities) adds only a few cents per kg.")


if __name__ == "__main__":
    main()
