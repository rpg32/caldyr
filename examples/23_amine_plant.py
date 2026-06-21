"""M16 demo: the FULL amine gas-sweetening plant (Hameed §15.3) — absorber +
regenerator + lean-amine recycle, plus the regenerator **reflux condenser**.

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

Part 2 then adds a **partial reflux condenser** to the regenerator top stage
(``condenser_T``): it condenses the overhead water back as internal reflux, so
the acid-gas product leaves *dried* (the §15.3 refinement). Drying the product
collapses the water carried off with the acid gas; closing the water loop to a
near-zero makeup additionally wants a reboiler in place of the open steam (a
reboiled-NS regenerator — endgame-stiff on this reactive system, tracked as
follow-up), so the validated recycle below uses open steam plus makeup.

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
T_COND = 320.0       # reflux-condenser temperature, K (~47 C)
Z_SOUR = {"methane": 0.9415, "CO2": 0.0413, "H2S": 0.0172}


def _stream(flows, T, P):
    F = sum(flows.values())
    return Stream(id="s", components=COMPS, T=T, P=P, molar_flow=F,
                  z={c: flows.get(c, 0.0) / F for c in COMPS})


def _flows(s):
    return {c: s.molar_flow * s.z[c] for c in COMPS}


def run_plant() -> dict:
    """Converge the absorber -> regenerator -> lean-recycle loop (open steam)."""
    pp = make_package("amine:DEA", COMPS)
    sour = {c: F_SOUR * Z_SOUR.get(c, 0.0) for c in COMPS}
    lean = {"DEA": 52.0, "water": W_CIRC, "CO2": 0.05, "H2S": 0.02, "methane": 0.0}

    print("Amine gas-sweetening plant (Hameed 2025 §15.3) — converging the loop")
    makeup = 0.0
    a_out = r_out = None
    for it in range(1, 16):
        absrb = Absorber("ABS", {"n_stages": 16, "murphree": {"CO2": 0.6}})
        a_out = absrb.solve({"gas_in": _stream(sour, 313.15, P_ABS),
                             "liquid_in": _stream(lean, 313.15, P_ABS)}, pp)
        rich = _flows(a_out["liquid_out"])

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

    sweet, acid = a_out["vapor_out"], r_out["vapor_out"]
    rich = _flows(a_out["liquid_out"])
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
    return {"pp": pp, "rich": rich}


def demo_reflux_condenser(pp, rich) -> None:
    """Part 2: the regenerator reflux condenser dries the acid-gas product.

    Run the same regenerator on the converged rich amine WITHOUT and WITH a
    partial reflux condenser, and contrast the overhead.
    """
    print("\nRegenerator reflux condenser (dry the acid-gas product)")
    feed = {"liquid_in": _stream(rich, 388.0, P_REG),
            "gas_in": _stream({"water": F_STEAM}, 396.0, P_REG)}

    wet = Absorber("REGwet", {"n_stages": 9, "method": "naphtali_sandholm",
                              "max_iter": 80})
    aw = wet.solve(feed, pp)["vapor_out"]

    dry = Absorber("REGdry", {"n_stages": 9, "method": "naphtali_sandholm",
                              "max_iter": 120, "condenser_T": T_COND})
    ad = dry.solve(feed, pp)["vapor_out"]
    qc = dry.design["condenser_duty"]

    w_wet = aw.molar_flow * aw.z["water"]
    w_dry = ad.molar_flow * ad.z["water"]
    print(f"  open steam, no condenser : acid gas {aw.molar_flow:6.1f} mol/s, "
          f"H2O {aw.z['water']:.2f}  (water carried off {w_wet:.1f} mol/s)")
    print(f"  + reflux condenser {T_COND:.0f} K : acid gas {ad.molar_flow:6.1f} mol/s, "
          f"H2O {ad.z['water']:.2f}  (water carried off {w_dry:.1f} mol/s)")
    print(f"  CO2 {ad.z['CO2']:.2f} / H2S {ad.z['H2S']:.2f} in the dried product; "
          f"condenser duty {qc / 1e6:.2f} MW")
    print(f"  -> water carried off with the acid gas cut {w_wet / w_dry:.0f}x "
          f"(the reflux returns it to the solvent).")
    print("  Note: with open steam the dried overhead leaves a steam-water "
          "surplus; a reboiler (no open steam) closes the loop to a small makeup.")


def main() -> None:
    state = run_plant()
    demo_reflux_condenser(state["pp"], state["rich"])


if __name__ == "__main__":
    main()
