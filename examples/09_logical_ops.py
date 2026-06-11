"""M8 demo: flowsheet-level logical ops — Set + Adjust — plus the balance report.

    FEED ─▶ HEAT ─▶ FL1 ─┬─▶ VAP1 (light product)
                         │
                 (liquid)▼
                        FL2 ─┬─▶ VAP2
                             └─▶ BOT (heavy product)

Two flash drums in series. The logical ops (stored on the flowsheet and
persisted in `.flow` under the "logical" key, like HYSYS SET/ADJUST blocks):

  * Set    — FL1.T is locked to HEAT.T_out (the first drum runs isothermal at
             whatever the heater delivers): FL1.T = 1.0*HEAT.T_out + 0.
  * Set    — FL2.T runs 15 K hotter than FL1.T (a reboil step, so the second
             drum actually makes vapor): FL2.T = 1.0*FL1.T + 15.
  * Adjust — vary HEAT.T_out within [325, 372] K (between the feed bubble
             point ~331 K and dew point ~378 K, so FL1 always keeps a liquid
             leg) until the heavy product BOT carries 4.0 mol/s, solved by a
             Brent root find around the full flowsheet solve.

Order (documented contract): Sets re-apply before *every* inner solve, so the
Adjust drags both flash temperatures along on each iteration. Afterwards,
`solver.balance_report` audits per-unit and overall mass/energy closure.

Run from the repo root:

    python examples/09_logical_ops.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.solver import balance_report  # noqa: E402
from caldyr.unitops import FlashDrum, Heater  # noqa: E402

P_ATM = 101325.0


def build() -> Flowsheet:
    fs = Flowsheet(
        components=[Component("n-pentane"), Component("n-octane")],
        property_package="thermo:PR",
    )
    fs.add(Heater("HEAT", {"T_out": 350.0}))
    fs.add(FlashDrum("FL1", {"P": P_ATM}))          # T comes from the Set
    fs.add(FlashDrum("FL2", {"P": P_ATM}))          # T comes from the other Set

    fs.feed("FEED", "HEAT:in1", T=300.0, P=P_ATM, molar_flow=10.0,
            z={"n-pentane": 0.5, "n-octane": 0.5})
    fs.connect("S1", "HEAT:out", "FL1:in1")
    fs.connect("VAP1", "FL1:vapor", None)
    fs.connect("L1", "FL1:liquid", "FL2:in1")
    fs.connect("VAP2", "FL2:vapor", None)
    fs.connect("BOT", "FL2:liquid", None)
    for u in ("HEAT", "FL1", "FL2"):
        fs.connect(f"Q_{u}", f"{u}:duty", None)

    fs.logical = [
        {"type": "set", "target": ["FL1", "T"], "source": ["HEAT", "T_out"]},
        {"type": "set", "target": ["FL2", "T"], "source": ["FL1", "T"], "offset": 15.0},
        {"type": "adjust", "vary": ["HEAT", "T_out"], "bounds": [325.0, 372.0],
         "spec": {"type": "molar_flow", "stream": "BOT"},
         "value": 4.0, "tolerance": 1e-7},
    ]
    return fs


def main() -> None:
    fs = build()
    report = fs.solve()

    print(f"Converged: {report.converged}")
    print("Solve / logical-op messages:")
    for msg in report.messages:
        print(f"  {msg}")

    print(f"\nAdjust result: HEAT.T_out = {fs.units['HEAT'].params['T_out']:.2f} K "
          f"-> FL1.T = {fs.units['FL1'].params['T']:.2f} K (Set), "
          f"FL2.T = {fs.units['FL2'].params['T']:.2f} K (Set)")

    head = f"{'stream':>6} {'T/K':>7} {'n/(mol/s)':>10} {'x_C5':>7}"
    print("\n" + head + "\n" + "-" * len(head))
    for sid in ("FEED", "VAP1", "VAP2", "BOT"):
        s = fs.streams[sid]
        print(f"{sid:>6} {s.T:>7.2f} {s.molar_flow:>10.4f} "
              f"{s.z.get('n-pentane', 0.0):>7.4f}")

    br = balance_report(fs)
    o = br["overall"]
    print("\nBalance report (worst offenders first):")
    print(f"  overall: mass {o['mass_in_kg_s']:.6f} -> {o['mass_out_kg_s']:.6f} kg/s "
          f"(rel {o['mass_rel']:.1e}); energy rel {o['energy_rel']:.1e} "
          f"(duties {o['duty_W'] / 1e3:+.2f} kW)")
    for u in br["units"]:
        print(f"  {u['unit_id']:>6}: mass rel {u['mass_rel']:.1e}, "
              f"energy rel {u['energy_rel']:.1e}, duty {u['duty_W'] / 1e3:+.2f} kW")
    for w in br["warnings"]:
        print(f"  WARNING: {w}")


if __name__ == "__main__":
    main()
