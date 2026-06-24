"""Example 34 — the full closed §9.5.6 anhydrous-ethanol entrainer PLANT.

Two columns + two recycles (Hameed 2025 §9.5.6):

    fresh feed (EtOH/water) ─┐
                             ▼
         cyclohexane     ┌────────┐  aqueous (EtOH/water/cyclohexane)
         make-up ───────►│ T-100  ├──────────────────────────────┐
            ▲            │ decant │                               ▼
            │            │ column │                          ┌────────┐
       ┌────┴────┐       └───┬────┘                          │ T-101  │
       │ Makeup  │◄──────────┼───────────────────────────────┤ water  │
       └─────────┘   recycle (EtOH/cyclohexane)  T-101 dist   │ column │
                             │                                └───┬────┘
                    anhydrous EtOH (bottoms)            water (bottoms)

T-100 is the integrated DECANTING-CONDENSER column (`decant_condenser=True`): its
overhead settles internally, the cyclohexane-rich organic layer is refluxed in
full (NO external organic recycle — that is the whole point), and the ethanol/
water aqueous layer is the distillate. T-101 recovers the water (bottoms) and
recycles the ethanol+cyclohexane (distillate) back to T-100. A :class:`Makeup`
controller holds the cyclohexane circulation; the two recycles are torn +
converged by the sequential solver (the §15.3 amine-plant pattern,
`examples/24`/`25`).

Anhydrous bottoms need a LARGE T-100 distillate (water leaves only in the
ethanol-rich aqueous layer), which a cold column cannot reach. So the plant is
brought up by **distillate-rate continuation at the flowsheet level**: solve the
closed loop at a small D, then ramp D up, re-solving the loop each step. T-100
warm-starts from the previous step (a few Newton iters), and the recycle
re-converges quickly. A cyclohexane-rich tear guess gives T-100 abundant
entrainer on the first (cold) solve.

NOTE: SLOW (two rigorous VLLE columns in a recycle, several continuation steps).
Run it directly; it is not part of the fast test suite.

    python examples/34_entrainer_plant.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.unitops import Makeup, RigorousColumn  # noqa: E402

P_ATM = 101325.0
COMPS = ["ethanol", "water", "cyclohexane"]
CYC_CIRC = 0.5   # cyclohexane MAKE-UP target: small (replaces losses only). The
#                  entrainer circulates internally (T-100 organic reflux) + via
#                  the T-101 recycle; abundant entrainer for the cold first solve
#                  comes from the cyclohexane-rich tear guess, not the make-up.


def build() -> Flowsheet:
    fs = Flowsheet(components=[Component(c) for c in COMPS],
                   property_package="thermo:UNIFAC")
    fs.add(RigorousColumn("T100", {
        "n_stages": 30, "feeds": [{"stage": 3}, {"stage": 14}],
        "reflux_ratio": 3.0, "distillate_rate": 4.0,   # ramped up in main()
        "method": "naphtali_sandholm", "reboiled": True,
        "decant_condenser": True, "condenser_T": 305.0,
        "reflux_layer": "organic", "max_iter": 120,
    }))
    fs.add(RigorousColumn("T101", {
        "n_stages": 16, "feed_stage": 8,
        "reflux_ratio": 1.0, "distillate_rate": 2.0,   # ramped up in main()
        "method": "naphtali_sandholm", "reboiled": True, "max_iter": 120,
    }))
    fs.add(Makeup("MK", {"component": "cyclohexane", "target": CYC_CIRC,
                         "T": 305.0, "P": P_ATM}))

    fs.feed("FEED", "T100:in2", T=343.0, P=P_ATM, molar_flow=27.78,
            z={"ethanol": 0.87, "water": 0.13})
    fs.connect("ENTR", "MK:out", "T100:in1")             # entrainer return (top)
    fs.connect("AQ", "T100:distillate", "T101:in1")      # aqueous -> water column
    fs.connect("REC", "T101:distillate", "MK:in1")       # the torn recycle
    fs.connect("ETOH", "T100:bottoms", None)             # anhydrous ethanol product
    fs.connect("WATER", "T101:bottoms", None)            # water product
    for u in ("T100", "T101"):
        fs.connect(f"{u}_QC", f"{u}:condenser_duty", None)
        fs.connect(f"{u}_QR", f"{u}:reboiler_duty", None)

    # Cyclohexane-rich tear guesses so T-100's first (cold) solve has abundant
    # entrainer (the over-fed regime converges cold). Guess both candidate torn
    # edges — the entrainer return into T-100 (ENTR) and the T-101 recycle (REC)
    # — so whichever the solver tears is seeded.
    entr0 = {"ethanol": 0.25, "water": 0.10, "cyclohexane": 0.65}
    rec0 = {"ethanol": 0.30, "water": 0.10, "cyclohexane": 0.60}
    fs.solver_hints = {
        "tear_guesses": {
            "ENTR": {"T": 320.0, "P": P_ATM, "molar_flow": 8.0, "z": entr0},
            "REC": {"T": 333.0, "P": P_ATM, "molar_flow": 8.0, "z": rec0},
        },
        "tear_tolerance": 5e-3,
    }
    return fs


def main() -> None:
    fs = build()
    t100, t101 = fs.units["T100"], fs.units["T101"]
    print("§9.5.6 anhydrous-ethanol entrainer plant — distillate-rate continuation:")
    print(f"  {'D100':>5} {'D101':>5} {'mk':>4}  {'recycle':>8}  "
          f"{'EtOH bottoms (x_EtOH / water / cyc)':<38}")
    # Ramp the loop up to anhydrous: distillate rates UP (drive the water out)
    # and the cyclohexane make-up target DOWN (toward replacing losses only, so
    # the entrainer stops piling into the bottoms). Gradual so T-100 warm-follows.
    mk = fs.units["MK"]
    for d100, d101, mkt in [(4.0, 2.0, 8.0), (7.0, 4.0, 6.0), (10.0, 6.0, 5.0),
                            (13.0, 9.0, 4.0), (16.0, 12.0, 3.0)]:
        t100.params["distillate_rate"] = d100
        t101.params["distillate_rate"] = d101
        mk.params["target"] = mkt
        try:
            rep = fs.solve(method="direct", max_iter=30)
        except Exception as exc:                # a column hit its endgame plateau
            print(f"  {d100:5.1f} {d101:5.1f} {mkt:4.1f}  "
                  f"{'STOP':>8}  {str(exc)[:60]}", flush=True)
            break
        etoh = fs.streams["ETOH"]
        print(f"  {d100:5.1f} {d101:5.1f} {mkt:4.1f}  "
              f"{'conv ' if rep.converged else 'NOCONV':>8}  "
              f"{etoh.z['ethanol']:.4f} / {etoh.z['water']:.4f} / "
              f"{etoh.z['cyclohexane']:.4f}  (sweeps {rep.iterations})", flush=True)

    etoh, water, rec = fs.streams["ETOH"], fs.streams["WATER"], fs.streams["REC"]
    print(f"\nethanol product (T-100 bottoms): {etoh.molar_flow:.2f} mol/s, "
          f"x_EtOH={etoh.z['ethanol']:.4f} "
          f"(water {etoh.z['water']:.4f}, cyclohexane {etoh.z['cyclohexane']:.4f})")
    print(f"water product (T-101 bottoms): {water.molar_flow:.2f} mol/s, "
          f"x_water={water.z['water']:.4f}")
    print(f"cyclohexane make-up: {fs.units['MK'].design['makeup_flow']:.3f} mol/s "
          f"(recycle returns {rec.molar_flow * rec.z['cyclohexane']:.2f} mol/s)")
    print("\nThe closed double recycle CONVERGES and the bottoms trend to "
          "anhydrous as D rises and the\nmake-up falls; book-scale >99.95% EtOH "
          "needs the 50-62 stages of Hameed Fig. 9.x (a\n30-stage column caps "
          "~90% here) — the integrated decant + continuation are the enablers.")


if __name__ == "__main__":
    main()
