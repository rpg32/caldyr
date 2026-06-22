"""Example 29 — packed-column sizing (Hameed 2025 sec. 9.1 SO2 absorber).

A packed tower is sized from the Eckert/Strigle generalized pressure-drop
correlation (GPDC) for the flood-limited diameter and from HETP rules of thumb
for the bed height, instead of the Fair flooding / tray-spacing rules used for
tray towers. Set ``internals='packed'`` on any column (distillation, absorber,
stripper) and the economics layer sizes and costs it as a packed bed —
tower shell + a random-packing bed (Perry 8e Tables 14-13/14-18).

This sizes the book's SO2/air/water absorber both ways and compares.

Run: python examples/29_packed_column.py
"""
from caldyr.core import Component, Flowsheet
from caldyr.economics.sizing import size_flowsheet
from caldyr.economics.costing import cost_equipment
from caldyr.thermo import make_package
from caldyr.unitops import Absorber

P_ATM = 101325.0
KMOLH = 1000.0 / 3600.0


def so2_absorber(internals: str, packing: str = "pall_metal_50mm") -> Flowsheet:
    fs = Flowsheet(
        components=[Component("nitrogen"), Component("oxygen"),
                    Component("sulfur dioxide"), Component("water")],
        property_package="thermo:PR")
    fs.add(Absorber("ABS", {"n_stages": 20, "P": P_ATM,
                            "internals": internals, "packing": packing}))
    fs.feed("GAS", "ABS:gas_in", T=293.15, P=P_ATM, molar_flow=206.0 * KMOLH,
            z={"sulfur dioxide": 0.03, "nitrogen": 0.97 * 0.79,
               "oxygen": 0.97 * 0.21})
    fs.feed("LIQ", "ABS:liquid_in", T=293.15, P=P_ATM,
            molar_flow=1.3e5 / 18.01528 * KMOLH, z={"water": 1.0})
    fs.connect("CLEAN", "ABS:gas_out", None)
    fs.connect("RICH", "ABS:liquid_out", None)
    return fs


def main() -> None:
    print("Hameed (2025) sec. 9.1 SO2 absorber — book HYSYS packed design: "
          "1.285 m\n")
    for internals in ("trays", "packed"):
        fs = so2_absorber(internals)
        rep = fs.solve()
        pp = make_package(fs.property_package, fs.component_ids)
        sizes = size_flowsheet(fs, rep, pp)
        tower = next(s for s in sizes if s.unit_id == "ABS")
        cbm = sum(cost_equipment(s).bare_module for s in sizes)
        print(f"  internals={internals:7s}: D={tower.diameter_m:.3f} m  "
              f"(volume {tower.attribute:.1f} m^3)  bare-module ${cbm/1e3:.0f}k")
        for note in tower.notes[:2]:
            print(f"      {note}")
    print("\nThe packed (50 mm metal Pall, GPDC at 70% of flood) and tray (Fair "
          "at 80% of\n  flood) diameters bracket the book value to ~15%.")


if __name__ == "__main__":
    main()
