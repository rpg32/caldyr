"""M16 demo: the §15.3 amine plant with a fully **reboiled + refluxed**
regenerator — a reflux condenser (dry acid-gas product) AND a reboiler duty (no
open steam), so the water loop closes and the make-up collapses to ~1 mol/s.

This is the closed-water-loop refinement of `examples/24` (which used open
stripping steam and a ~84 mol/s make-up). The regenerator here has:
  * ``condenser_T`` — a partial reflux condenser that dries the overhead;
  * ``reboiler_duty`` — the reboiler heat that boils up the stripping vapour
    internally (no ``gas_in``: no open steam, so no net water is added);
so the lean amine keeps ~all its water and the only water leaving the system is
the small amount in the dried acid-gas product.

Pinning both end temperatures (condenser T + reboiler T) is FD-Jacobian-stiff, so
the Absorber uses ``reboiler_duty`` and an internal warm-start continuation (a
robust open-steam proxy seeds the real internally-boiled solve) — it converges in
a handful of Newton steps. The lean-amine recycle is torn by the sequential
solver and the circulating water held by the in-loop :class:`Makeup` controller,
exactly as in `examples/24`.

    python examples/25_amine_plant_reboiled.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.unitops import Absorber, Heater, Makeup  # noqa: E402

COMPS = ["DEA", "water", "CO2", "H2S", "methane"]
P_ABS, P_REG = 6.9e6, 1.8e5
W_CIRC = 800.0
Q_REB = 9.0e6                 # regenerator reboiler duty, W
Z_SOUR = {"methane": 0.9415, "CO2": 0.0413, "H2S": 0.0172}


def build() -> Flowsheet:
    fs = Flowsheet(components=[Component(c) for c in COMPS],
                   property_package="amine:DEA")
    fs.add(Absorber("ABS", {"n_stages": 16, "murphree": {"CO2": 0.6}}))
    fs.add(Heater("HEAT", {"T_out": 388.0, "dP": P_ABS - P_REG}))      # heat + let down
    # Reboiled + refluxed regenerator: dry product (condenser) + internal boilup
    # (reboiler duty), NO open steam -> the water loop closes.
    fs.add(Absorber("REG", {"n_stages": 10, "method": "naphtali_sandholm",
                            "max_iter": 100, "condenser_T": 320.0,
                            "reboiler_duty": Q_REB}))
    fs.add(Heater("COOL", {"T_out": 313.15, "dP": P_REG - P_ABS}))     # cool + pump up
    fs.add(Makeup("MK", {"component": "water", "target": W_CIRC,
                         "T": 313.15, "P": P_ABS}))

    fs.feed("SOUR", "ABS:gas_in", T=313.15, P=P_ABS, molar_flow=350.0,
            z={c: Z_SOUR.get(c, 0.0) for c in COMPS})
    fs.connect("RICH", "ABS:liquid_out", "HEAT:in1")
    fs.connect("RICHHOT", "HEAT:out", "REG:liquid_in")
    fs.connect("LEAN", "REG:liquid_out", "COOL:in1")
    fs.connect("LEANCOOL", "COOL:out", "MK:in1")
    fs.connect("LEANRECY", "MK:out", "ABS:liquid_in")                 # the torn recycle
    fs.connect("SWEET", "ABS:vapor_out", None)
    fs.connect("ACID", "REG:vapor_out", None)
    fs.connect("HQ", "HEAT:duty", None)
    fs.connect("CQ", "COOL:duty", None)

    lean0 = {"DEA": 52.0, "water": W_CIRC, "CO2": 0.05, "H2S": 0.02}
    tot = sum(lean0.values())
    fs.solver_hints = {
        "tear_guesses": {"LEANRECY": {"T": 313.15, "P": P_ABS, "molar_flow": tot,
                                      "z": {c: lean0.get(c, 0.0) / tot for c in COMPS}}},
        "tear_tolerance": 3e-3,
    }
    return fs


def main() -> None:
    fs = build()
    print("§15.3 amine plant — reboiled + refluxed regenerator (closed water loop)")
    rep = fs.solve(method="direct", max_iter=40)
    print(f"  recycle {('CONVERGED' if rep.converged else 'did NOT converge')} "
          f"in {rep.iterations} sweeps (residual {rep.residual:.2e})")

    fl = lambda s, c: s.molar_flow * s.z.get(c, 0.0)                   # noqa: E731
    lr, sweet, acid = fs.streams["LEANRECY"], fs.streams["SWEET"], fs.streams["ACID"]
    d = fs.units["REG"].design
    print(f"  make-up water {fs.units['MK'].design['makeup_flow']:.1f} mol/s "
          f"(vs ~84 for the open-steam stripper) holds circulating water "
          f"{fl(lr, 'water'):.0f} mol/s")
    print(f"  regenerator: reboiler {d['reboiler_duty'] / 1e6:.1f} MW, "
          f"condenser {d['condenser_duty'] / 1e6:.1f} MW")
    print(f"  lean loading  CO2 {lr.z['CO2'] / lr.z['DEA']:.4f}  "
          f"H2S {lr.z['H2S'] / lr.z['DEA']:.4f}")
    rem = 1.0 - fl(sweet, "CO2") / (350.0 * Z_SOUR["CO2"])
    print(f"  sweet gas {sweet.molar_flow:.1f} mol/s, CO2 {1e6 * sweet.z['CO2']:.0f} "
          f"ppm ({100 * rem:.1f}% CO2 removed)")
    print(f"  acid gas {acid.molar_flow:.1f} mol/s, H2O {acid.z['water']:.2f} "
          f"(dried; CO2 {acid.z['CO2']:.2f} / H2S {acid.z['H2S']:.2f})")


if __name__ == "__main__":
    main()
