"""Example 32 — fired-heater radiant/convective design split (Hameed 2025 §4.3).

The book's §4.3.4 problem heats a 1250 kmol/h light-hydrocarbon (C1-C5) feed from
50 to 300 °C in a direct-fired heater burning methane at 100 % excess air and
60 % fired efficiency, and asks for the fuel and air rates and the flue-gas
temperature (book answer: ~62 kmol/h fuel, Fig. 4.19).

Caldyr's `FiredHeater` solves the process side as before; setting a fuel/excess-
air spec additionally runs the firebox combustion and the **radiant/convective
design split** (Hameed eqs. 4.8-4.9; API/Lobo-Evans furnace heat balance):
the fuel and air molar flows, the flue-gas composition, the adiabatic flame /
bridgewall / stack temperatures, the radiant-vs-convective duty division, and the
radiant and convective tube areas — all on `unit.design`.

Run: python examples/32_fired_heater_design.py
"""
from caldyr.core import Component, Flowsheet
from caldyr.unitops import FiredHeater

COMPS = ["methane", "ethane", "propane", "isobutane",
         "n-butane", "isopentane", "n-pentane"]
FLOWS = {"methane": 300, "ethane": 250, "propane": 200, "isobutane": 150,
         "n-butane": 150, "isopentane": 100, "n-pentane": 100}   # kmol/h


def main() -> None:
    tot = sum(FLOWS.values())                       # 1250 kmol/h
    z = {c: FLOWS[c] / tot for c in COMPS}
    fs = Flowsheet(components=[Component(c) for c in COMPS],
                   property_package="thermo:PR")
    fs.add(FiredHeater("FH", {
        "T_out": 300 + 273.15,
        "efficiency": 0.60,
        "fuel": "methane",          # triggers the radiant/convective design
        "excess_air": 1.0,          # 100 % excess air
        "fuel_T": 40 + 273.15,
        "air_T": 40 + 273.15,
    }))
    fs.feed("FEED", "FH:in1", T=50 + 273.15, P=5 * 101325.0,
            molar_flow=tot * 1000 / 3600.0, z=z)
    fs.connect("OUT", "FH:out", None)
    fs.connect("Q", "FH:duty", None)
    fs.solve()

    d = fs.units["FH"].design
    c, f = d["combustion"], d["firing"]
    print("Fired heater — Hameed §4.3 (1250 kmol/h C1-C5, 50 -> 300 °C):\n")
    print(f"  process (absorbed) duty   {d['process_duty'] / 1e6:7.3f} MW")
    print(f"  fired duty (LHV, Q/eta)   {f['fired_duty'] / 1e6:7.3f} MW")
    print(f"  fuel (methane)            {c['fuel_flow'] * 3.6:7.2f} kmol/h"
          f"   (book ~62 kmol/h; PR duty runs ~6 % high)")
    print(f"  combustion air            {c['air_flow'] * 3.6:7.1f} kmol/h")
    print(f"  flue gas                  {c['flue_flow'] * 3.6:7.1f} kmol/h"
          f"  {{{', '.join(f'{k.split()[0]}:{v:.2f}' for k, v in c['flue_composition'].items())}}}")
    print()
    print(f"  adiabatic flame T         {f['flame_temperature'] - 273.15:7.0f} °C")
    print(f"  bridgewall T (firebox)    {f['bridgewall_temperature'] - 273.15:7.0f} °C")
    print(f"  stack T                   {f['stack_temperature'] - 273.15:7.0f} °C")
    print()
    print(f"  radiant duty              {f['radiant_duty'] / 1e6:7.3f} MW"
          f"   ({f['radiant_fraction']:.0%} of absorbed)   area {f['radiant_area']:6.1f} m²")
    print(f"  convective duty           {f['convective_duty'] / 1e6:7.3f} MW"
          f"   (LMTD {f['convective_lmtd']:.0f} K)        area {f['convective_area']:6.1f} m²")
    print(f"  casing/stack loss         {(f['casing_loss'] + f['stack_loss']) / 1e6:7.3f} MW"
          f"   (gross efficiency {f['efficiency_gross']:.1%})")
    print("\n  Energy balance: Q_available = Q_absorbed + casing loss + stack loss")
    print(f"    {f['heat_available'] / 1e6:.3f} MW = {d['process_duty'] / 1e6:.3f}"
          f" + {f['casing_loss'] / 1e6:.3f} + {f['stack_loss'] / 1e6:.3f} MW")


if __name__ == "__main__":
    main()
