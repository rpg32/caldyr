"""M16 demo: the §15.3 amine gas-sweetening plant as a **native Flowsheet** —
the recycle solver tears and converges the lean-amine loop automatically, and
an in-loop make-up controller holds the circulating water.

This is the native-flowsheet form of `examples/23` (which converged the same
loop with a hand-rolled script-level tear). Here every block is a unit op on a
`Flowsheet`, and `fs.solve()` does the tearing:

    SOUR ->[ABS]-> SWEET (product)            STEAM ->[REG]
            |  ^                                      |  ^
       rich |  | lean (cooled, water made up)   acid v  | rich (heated)
            v  |                                  ACID    |
          [HEAT]---------------> RICHHOT ------------> REG ]
          [COOL]<--------------- LEAN <----------------'
            |
          [MK]  (Makeup: top water back to the circulating target)
            '--> LEANRECY --> ABS:liquid_in   (the torn recycle)

**Why a Makeup controller, not an outer Adjust.** An open-steam stripper loses
water overhead, so the lean-amine inventory is *unstable* under a fixed make-up:
the circulating water runs away (no stable fixed point), and an outer Adjust that
varies a make-up Source is fragile — at an off-target make-up the columns are
driven to states they cannot even solve. The robust realization is the
:class:`~caldyr.unitops.Makeup` controller, which tops the circulating water back
to its target *analytically each sweep* (the make-up rate it computes IS what an
Adjust would search for) — so the recycle's water mode is pinned at its fixed
point and the loop converges in a few direct-substitution sweeps.

Direct substitution (not Wegstein) is used: the inner column solvers inject a
little noise into the tear, which destabilizes Wegstein's slope estimate; with
the water pinned, plain substitution converges cleanly. The tear tolerance is
~3e-3 because the *lean* loading is ~3e-4, so the trace acid-gas flows have a
large relative sensitivity while the macroscopic balances are converged.

    python examples/24_amine_plant_flowsheet.py
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
Z_SOUR = {"methane": 0.9415, "CO2": 0.0413, "H2S": 0.0172}


def build() -> Flowsheet:
    fs = Flowsheet(components=[Component(c) for c in COMPS],
                   property_package="amine:DEA")
    fs.add(Absorber("ABS", {"n_stages": 16, "murphree": {"CO2": 0.6}}))
    fs.add(Heater("HEAT", {"T_out": 388.0, "dP": P_ABS - P_REG}))   # heat + let down
    fs.add(Absorber("REG", {"n_stages": 8, "method": "naphtali_sandholm",
                            "max_iter": 80}))
    fs.add(Heater("COOL", {"T_out": 313.15, "dP": P_REG - P_ABS}))  # cool + pump up
    fs.add(Makeup("MK", {"component": "water", "target": W_CIRC,
                         "T": 313.15, "P": P_ABS}))

    fs.feed("SOUR", "ABS:gas_in", T=313.15, P=P_ABS, molar_flow=350.0,
            z={c: Z_SOUR.get(c, 0.0) for c in COMPS})
    fs.feed("STEAM", "REG:gas_in", T=396.0, P=P_REG, molar_flow=140.0,
            z={"water": 1.0})
    fs.connect("RICH", "ABS:liquid_out", "HEAT:in1")
    fs.connect("RICHHOT", "HEAT:out", "REG:liquid_in")
    fs.connect("LEAN", "REG:liquid_out", "COOL:in1")
    fs.connect("LEANCOOL", "COOL:out", "MK:in1")
    fs.connect("LEANRECY", "MK:out", "ABS:liquid_in")              # the torn recycle
    fs.connect("SWEET", "ABS:vapor_out", None)
    fs.connect("ACID", "REG:vapor_out", None)
    fs.connect("HQ", "HEAT:duty", None)
    fs.connect("CQ", "COOL:duty", None)

    # Seed the lean-amine recycle so the absorber's first sweep has solvent.
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
    print("§15.3 amine plant as a native Flowsheet — solving the recycle")
    rep = fs.solve(method="direct", max_iter=40)
    print(f"  recycle {('CONVERGED' if rep.converged else 'did NOT converge')} "
          f"in {rep.iterations} sweeps (residual {rep.residual:.2e}); "
          f"torn stream {rep.tear_streams}")

    fl = lambda s, c: s.molar_flow * s.z.get(c, 0.0)                   # noqa: E731
    lr, rich = fs.streams["LEANRECY"], fs.streams["RICH"]
    sweet, acid = fs.streams["SWEET"], fs.streams["ACID"]
    print(f"  make-up water {fs.units['MK'].design['makeup_flow']:.1f} mol/s "
          f"holds circulating water {fl(lr, 'water'):.0f} mol/s, DEA {fl(lr, 'DEA'):.1f}")
    print(f"  rich loading  CO2 {rich.z['CO2'] / rich.z['DEA']:.3f}  "
          f"H2S {rich.z['H2S'] / rich.z['DEA']:.3f}")
    print(f"  lean loading  CO2 {lr.z['CO2'] / lr.z['DEA']:.4f}  "
          f"H2S {lr.z['H2S'] / lr.z['DEA']:.4f}")
    rem = 1.0 - fl(sweet, "CO2") / (350.0 * Z_SOUR["CO2"])
    print(f"  sweet gas {sweet.molar_flow:.1f} mol/s, CO2 {1e6 * sweet.z['CO2']:.0f} "
          f"ppm, H2S {1e6 * sweet.z['H2S']:.0f} ppm ({100 * rem:.1f}% CO2 removed)")
    print(f"  acid gas {acid.molar_flow:.1f} mol/s "
          f"(CO2 {acid.z['CO2']:.2f} / H2S {acid.z['H2S']:.2f} / H2O {acid.z['water']:.2f})")

    # Overall plant component balance: feeds + make-up water == products.
    makeup = fs.units["MK"].design["makeup_flow"]
    sour, steam = fs.streams["SOUR"], fs.streams["STEAM"]
    worst = 0.0
    for c in COMPS:
        cin = fl(sour, c) + fl(steam, c) + (makeup if c == "water" else 0.0)
        cout = fl(sweet, c) + fl(acid, c)
        worst = max(worst, abs(cin - cout))
    print(f"  overall component balance (feeds + make-up = products) closes to "
          f"{worst:.2e} mol/s")


if __name__ == "__main__":
    main()
