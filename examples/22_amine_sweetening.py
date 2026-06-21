"""M16 demo: amine gas sweetening — a DEA absorber over the reactive acid-gas package.

Scrubs CO2 and H2S out of a sour natural gas with aqueous DEA, after Hameed,
"Chemical Process Simulations using Aspen HYSYS" (Wiley 2025), sec. 15.3 "Gas
Sweetening". The chemistry that makes this work — acid gas reacting with the amine
in the liquid to form ionic species, so only the small free-molecular fraction has
a vapour pressure — is captured by the modified Kent-Eisenberg property package
(`property_package = "amine:DEA"`, the open analogue of HYSYS's "Acid Gas -
Chemical Solvents"). The acid-gas K-value collapses at low loading (avid
absorption) and rises as the solvent saturates.

The column is a plain equilibrium-stage absorber (no condenser/reboiler), solved
by the sum-rates (Burningham-Otto) method. The book gives no numeric output slate,
so the demo reports the physics: acid-gas removal, the rich loading, and the
temperature bulge from the heat of absorption.

    python examples/22_amine_sweetening.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.core import Stream  # noqa: E402
from caldyr.thermo import make_package  # noqa: E402
from caldyr.unitops.absorber import Absorber  # noqa: E402

COMPS = ["DEA", "water", "CO2", "H2S", "methane"]
P = 6.9e6        # ~1000 psia contactor pressure


def main() -> None:
    pp = make_package("amine:DEA", COMPS)

    # Sour natural gas: ~4.1% CO2, ~1.7% H2S, balance methane (book §15.3
    # proportions); ~350 mol/s ~ 25 MMSCFD.
    z_g = {"methane": 0.9415, "CO2": 0.0413, "H2S": 0.0172}
    sour = Stream(id="sour_gas", components=COMPS, T=313.15, P=P, molar_flow=350.0,
                  z={c: z_g.get(c, 0.0) for c in COMPS})

    # Lean solvent: 28 wt% (~2.7 M) aqueous DEA with a small residual loading.
    z_l = {"DEA": 0.062, "water": 0.9355, "CO2": 0.0017, "H2S": 0.0008}
    lean = Stream(id="lean_amine", components=COMPS, T=313.15, P=P, molar_flow=850.0,
                  z={c: z_l.get(c, 0.0) for c in COMPS})

    column = Absorber("contactor", {"n_stages": 20})
    out = column.solve({"gas_in": sour, "liquid_in": lean}, pp)
    d = column.design
    sweet = out["vapor_out"]

    print("DEA gas-sweetening absorber (Hameed 2025 §15.3)")
    print(f"  converged in {d['iterations']} sum-rates iterations "
          f"(energy residual {d['energy_residual_rel']:.1e})")
    print(f"  acid-gas removal:  CO2 {100 * d['absorbed']['CO2']:.1f}%   "
          f"H2S {100 * d['absorbed']['H2S']:.1f}%   "
          f"(methane slip {100 * d['absorbed']['methane']:.2f}%)")
    print(f"  sweet gas: {sweet.molar_flow:.1f} mol/s, "
          f"CO2 {1e6 * sweet.z['CO2']:.0f} ppm, H2S {1e6 * sweet.z['H2S']:.0f} ppm")
    xb = d["x_bottom"]
    print(f"  rich loading:  CO2 {xb['CO2'] / xb['DEA']:.3f}   "
          f"H2S {xb['H2S'] / xb['DEA']:.3f}  mol/mol amine")
    print(f"  temperature bulge: {d['T_top'] - 273.15:.1f} -> "
          f"{d['T_bottom'] - 273.15:.1f} C (heat of absorption)")

    _mdea_selective_demo()


def _mdea_selective_demo() -> None:
    """MDEA selectively removes H2S over CO2 — but the selectivity is *kinetic*
    (CO2 reacts slowly with the tertiary amine), so it only appears once a
    Murphree stage efficiency limits the CO2 mass transfer (book §15.3 uses
    E_CO2 ~ 0.15, E_H2S ~ 0.8). At full equilibrium an MDEA stage absorbs both."""
    comps = ["MDEA", "water", "CO2", "H2S", "methane"]
    pp = make_package("amine:MDEA", comps)
    z_g = {"methane": 0.9415, "CO2": 0.0413, "H2S": 0.0172}
    z_l = {"MDEA": 0.062, "water": 0.9355, "CO2": 0.0017, "H2S": 0.0008}

    def removal(murphree):
        gas = Stream(id="g", components=comps, T=313.15, P=P, molar_flow=350.0,
                     z={c: z_g.get(c, 0.0) for c in comps})
        liq = Stream(id="l", components=comps, T=313.15, P=P, molar_flow=600.0,
                     z={c: z_l.get(c, 0.0) for c in comps})
        params = {"n_stages": 16}
        if murphree:
            params["murphree"] = murphree
        col = Absorber("mdea", params)
        col.solve({"gas_in": gas, "liquid_in": liq}, pp)
        return col.design["absorbed"]

    eq = removal(None)
    kin = removal({"CO2": 0.15, "H2S": 0.8})
    print("\nMDEA selectivity (the role of Murphree efficiency)")
    print(f"  equilibrium stages:  CO2 {100 * eq['CO2']:.0f}%  H2S {100 * eq['H2S']:.0f}%"
          f"   (no selectivity)")
    print(f"  E_CO2=0.15/E_H2S=0.8: CO2 {100 * kin['CO2']:.0f}%  H2S {100 * kin['H2S']:.0f}%"
          f"   (H2S-selective, as in practice)")


if __name__ == "__main__":
    main()
