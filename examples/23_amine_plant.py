"""M16 demo: the FULL amine gas-sweetening plant (Hameed §15.3) — absorber +
regenerator + lean-amine recycle.

Builds on examples/22 (the standalone DEA absorber) to close the loop:

    sour gas ->[ABSORBER]-> sweet gas
                  |  ^
            rich  |  | lean amine (recycle, cooled)
                  v  |
        (heat to ~115 C) ->[REGENERATOR]-> acid gas (overhead)
                                |
                          lean amine (hot, bottoms)

The **absorber** is the sum-rates equilibrium column with a Murphree efficiency
on CO2 (kinetic limitation). The **regenerator** is a steam stripper at ~1.8 bar
solved by the Naphtali-Sandholm simultaneous-correction method — sum-rates
limit-cycles on reactive desorption and the bubble-point method hits a degenerate
energy balance, so NS is required here. The lean amine is recycled; a small
water makeup holds the circulating water (open-steam stripping carries more water
overhead than the steam supplies, exactly as a real unit's makeup compensates).

The recycle tear (the lean amine) is converged by damped direct substitution.

    python examples/23_amine_plant.py
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
P_ABS = 6.9e6        # contactor ~1000 psia
P_REG = 1.8e5        # regenerator ~1.8 bar
F_SOUR = 350.0       # sour gas, mol/s (~25 MMSCFD)
F_STEAM = 140.0      # stripping steam, mol/s
W_CIRC = 800.0       # circulating water held by makeup, mol/s
Z_SOUR = {"methane": 0.9415, "CO2": 0.0413, "H2S": 0.0172}


def _stream(flows, T, P):
    F = sum(flows.values())
    return Stream(id="s", components=COMPS, T=T, P=P, molar_flow=F,
                  z={c: flows.get(c, 0.0) / F for c in COMPS})


def _flows(s):
    return {c: s.molar_flow * s.z[c] for c in COMPS}


def main() -> None:
    pp = make_package("amine:DEA", COMPS)
    sour = {c: F_SOUR * Z_SOUR.get(c, 0.0) for c in COMPS}
    lean = {"DEA": 52.0, "water": W_CIRC, "CO2": 0.05, "H2S": 0.02, "methane": 0.0}

    print("Amine gas-sweetening plant (Hameed 2025 §15.3) — converging the loop")
    makeup = 0.0
    for it in range(1, 16):
        absrb = Absorber("ABS", {"n_stages": 16, "murphree": {"CO2": 0.6}})
        a_out = absrb.solve({"gas_in": _stream(sour, 313.15, P_ABS),
                             "liquid_in": _stream(lean, 313.15, P_ABS)}, pp)
        rich = _flows(a_out["liquid_out"])
        sweet = a_out["vapor_out"]

        regen = Absorber("REG", {"n_stages": 8, "method": "naphtali_sandholm",
                                 "max_iter": 80})
        r_out = regen.solve({"liquid_in": _stream(rich, 388.0, P_REG),
                             "gas_in": _stream({"water": F_STEAM}, 396.0, P_REG)},
                            pp)
        new = _flows(r_out["liquid_out"])
        makeup = W_CIRC - new["water"]            # water makeup closes the balance
        new["water"] = W_CIRC
        err = sum(abs(new[c] - lean[c]) for c in COMPS) / sum(lean.values())
        lean = {c: 0.5 * lean[c] + 0.5 * new[c] for c in COMPS}
        print(f"  iter {it:2d}: tear residual {err:.2e}")
        if err < 1e-4:
            break

    acid = r_out["vapor_out"]
    co2_rem = 1.0 - sweet.z["CO2"] * sweet.molar_flow / sour["CO2"]
    h2s_rem = 1.0 - sweet.z["H2S"] * sweet.molar_flow / sour["H2S"]
    print("\nConverged plant")
    print(f"  circulating DEA {lean['DEA']:.1f} mol/s, water {lean['water']:.0f} "
          f"mol/s; makeup water {makeup:.1f} mol/s")
    print(f"  rich loading  CO2 {rich['CO2'] / rich['DEA']:.3f}  "
          f"H2S {rich['H2S'] / rich['DEA']:.3f}")
    print(f"  lean loading  CO2 {lean['CO2'] / lean['DEA']:.4f}  "
          f"H2S {lean['H2S'] / lean['DEA']:.4f}")
    print(f"  acid-gas removal  CO2 {100 * co2_rem:.1f}%  H2S {100 * h2s_rem:.1f}%")
    print(f"  sweet gas {sweet.molar_flow:.1f} mol/s, "
          f"CO2 {1e6 * sweet.z['CO2']:.0f} ppm, H2S {1e6 * sweet.z['H2S']:.0f} ppm")
    print(f"  acid-gas product {acid.molar_flow:.1f} mol/s "
          f"(CO2 {acid.z['CO2']:.2f} / H2S {acid.z['H2S']:.2f} / "
          f"H2O {acid.z['water']:.2f})")


if __name__ == "__main__":
    main()
