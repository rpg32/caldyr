"""Example 31 — gas saturator (Hameed 2025 sec. 10.4 Stream Saturator).

HYSYS's Stream Saturator humidifies the acid-gas and combustion-air feeds of a
sulfur-recovery unit — "used to quickly calculate saturation of acid gas and air
streams". The `Saturator` unit op loads a gas with a condensable component
(water by default) up to its saturation partial pressure at the operating
temperature and pressure (or a chosen relative humidity); the latent heat shows
up honestly on the duty port.

This saturates a dry nitrogen stream with n-hexane (the cubic EOS predicts the
hexane vapour pressure within a few percent — water VLE under PR is poor, the
documented caveat) at 300 K and 1 atm, then shows the relative-humidity knob and
an excess-liquid case.

Run: python examples/31_saturator.py
"""
from caldyr.core import Component, Flowsheet
from caldyr.economics.costing import cost_equipment
from caldyr.economics.sizing import size_flowsheet
from caldyr.thermo import make_package
from caldyr.unitops import Saturator

P_ATM = 101325.0
COMPS = ["nitrogen", "hexane"]


def run(rh, water=None):
    fs = Flowsheet(components=[Component(c) for c in COMPS],
                   property_package="thermo:PR")
    fs.add(Saturator("SAT", {"saturant": "hexane", "relative_humidity": rh}))
    fs.feed("GAS", "SAT:gas_in", T=300.0, P=P_ATM, molar_flow=1.0,
            z={"nitrogen": 1.0, "hexane": 0.0})
    if water is not None:
        fs.feed("LIQ", "SAT:water_in", T=300.0, P=P_ATM, molar_flow=water,
                z={"nitrogen": 0.0, "hexane": 1.0})
    fs.connect("WETGAS", "SAT:gas_out", None)
    fs.connect("DRAIN", "SAT:liquid_out", None)
    fs.connect("Q", "SAT:duty", None)
    rep = fs.solve()
    d = fs.units["SAT"].design
    g = fs.streams["WETGAS"]
    pp = make_package(fs.property_package, fs.component_ids)
    cbm = sum(cost_equipment(s).bare_module for s in size_flowsheet(fs, rep, pp))
    drain = fs.streams["DRAIN"].molar_flow
    print(f"  RH target {rh:.0%}: y_hexane={g.z['hexane']:.4f} "
          f"(p={g.z['hexane'] * g.P / 1e3:.1f} kPa), RH achieved "
          f"{d['relative_humidity_achieved']:.0%}, hexane added "
          f"{d['saturant_added']:.4f} mol/s, drain {drain:.3f} mol/s, "
          f"duty {d['duty'] / 1e3:.2f} kW, vessel ${cbm / 1e3:.0f}k")


def main() -> None:
    from thermo import Chemical
    psat = Chemical("hexane", T=300.0, P=P_ATM).Psat
    print("Saturate dry N2 with n-hexane at 300 K, 1 atm.")
    print(f"  (n-hexane vapour pressure at 300 K ~ {psat / 1e3:.1f} kPa "
          f"-> saturation y ~ {psat / P_ATM:.4f})\n")
    run(1.0)                       # fully saturated, auto water
    run(0.5)                       # half relative humidity
    print()
    run(1.0, water=0.40)           # excess hexane supply -> surplus drains
    print("\n  The saturated hexane partial pressure tracks its vapour pressure; "
          "relative\n  humidity scales it linearly; surplus liquid drains and the "
          "latent heat\n  to hold temperature is reported on the duty port.")


if __name__ == "__main__":
    main()
