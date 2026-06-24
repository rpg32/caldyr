"""Example 33 — anhydrous ethanol via an integrated DECANTING-CONDENSER column,
at BOOK SCALE (62 stages).

Hameed (2025) §9.5.6 dehydrates ethanol past its water azeotrope with a
**cyclohexane** entrainer. The trick is a heteroazeotropic column whose overhead
is the ternary (ethanol/water/cyclohexane) heteroazeotrope: it condenses and
settles into two liquids — a cyclohexane-rich ORGANIC layer and an ethanol/water
AQUEOUS layer.

The decisive modelling point (proven the hard way — see ``PROGRESS.md`` P6):

  * An **external** decanter feeds the entrainer back to the column as a stream,
    so the column overhead must carry the *entire* entrainer circulation out the
    top → a very large distillate the Naphtali-Sandholm solver cannot converge.
  * An **integrated decanting condenser** settles the overhead INSIDE the unit and
    refluxes the organic layer internally, so the cyclohexane never leaves — only
    the (small, net) aqueous layer is drawn. That keeps the entrainer circulating
    and makes the column tractable.

This is the ``RigorousColumn`` ``decant_condenser=True`` mode: stage 0 is an
ordinary tray whose overhead vapour is flashed at ``condenser_T`` via the
predictive UNIFAC VLLE package (``thermo:UNIFAC``), the organic layer is refluxed
in full, and the aqueous layer is the distillate (spec = ``distillate_rate``).

Reaching anhydrous bottoms needs (a) a long stripping section to strip the
cyclohexane (bp 81 ≈ ethanol 78 °C) out of the bottoms — Hameed uses ~62 stages —
and (b) a large aqueous draw, since water leaves ONLY in the ethanol-rich aqueous
layer. Two continuations get there:

  * **Stage-count continuation.** A cold 62-stage decant solve is intractable, so
    we cold-solve a cheap 30-stage column (~1 min) and ``warm_start_from`` it —
    the converged profile, interpolated onto 62 stages, lands the long column in
    its basin (it then re-solves in a few warm Newton iterations).
  * **Distillate-rate continuation.** From the seeded 62-stage column we ramp the
    aqueous draw D up (trimming the cyclohexane solvent feed down with it), each
    step warm-starting from the previous. By D≈18 the bottoms is ~99 % ethanol,
    water essentially gone (~1e-4), with the residual cyclohexane the last
    impurity (driven out by the inventory control closed in example 34).

Run: python examples/33_entrainer_decant_column.py
"""
from caldyr.core import Component, Flowsheet
from caldyr.unitops import RigorousColumn

P_ATM = 101325.0
COMPS = ["ethanol", "water", "cyclohexane"]


def _solve(col: RigorousColumn, D: float, solv: float):
    col.params["distillate_rate"] = D
    fs = Flowsheet(components=[Component(c) for c in COMPS],
                   property_package="thermo:UNIFAC")
    fs.add(col)
    fs.feed("SOLV", "T100:in1", T=298.15, P=P_ATM, molar_flow=solv,
            z={"cyclohexane": 1.0})
    fs.feed("FEED", "T100:in2", T=343.0, P=P_ATM, molar_flow=27.8,
            z={"ethanol": 0.87, "water": 0.13})
    fs.connect("DIST", "T100:distillate", None)
    fs.connect("BOT", "T100:bottoms", None)
    fs.connect("QC", "T100:condenser_duty", None)
    fs.connect("QR", "T100:reboiler_duty", None)
    fs.solve()
    return fs


def main() -> None:
    # -- step 1: cold 30-stage solve (the cheap basin) --------------------------
    print("Cold-solving a 30-stage proxy column (stage-count continuation seed)...",
          flush=True)
    col30 = RigorousColumn("T100", {
        "n_stages": 30, "feeds": [{"stage": 2}, {"stage": 10}],
        "reflux_ratio": 3.0, "distillate_rate": 4.2,
        "method": "naphtali_sandholm", "reboiled": True,
        "decant_condenser": True, "condenser_T": 305.0,
        "reflux_layer": "organic", "max_iter": 120,
    })
    _solve(col30, 4.2, 8.0)

    # -- step 2: seed a 62-stage column from it, then ramp D up -----------------
    col = RigorousColumn("T100", {
        "n_stages": 62, "feeds": [{"stage": 2}, {"stage": 20}],
        "reflux_ratio": 3.0, "distillate_rate": 4.2,
        "method": "naphtali_sandholm", "reboiled": True,
        "decant_condenser": True, "condenser_T": 305.0,
        "reflux_layer": "organic", "max_iter": 150,
    })
    col.warm_start_from(col30)

    print("\n62-stage integrated decanting-condenser column "
          "(ethanol/water/cyclohexane):")
    print(f"  {'D':>5} {'solv':>5} {'iters':>6} {'Eres':>9}   "
          f"{'bottoms ethanol / water / cyclohexane':<40}")
    for D, solv in [(4.2, 8.0), (6.0, 7.0), (8.0, 6.0), (10.0, 5.0), (12.0, 4.5),
                    (14.0, 4.0), (16.0, 3.6), (18.0, 3.0)]:
        fs = _solve(col, D, solv)
        d = col.design
        bot = fs.streams["BOT"]
        print(f"  {D:5.1f} {solv:5.1f} {d['iterations']:6d} "
              f"{d['energy_residual_rel']:9.1e}   "
              f"{bot.z['ethanol']:.5f} / {bot.z['water']:.5f} / "
              f"{bot.z['cyclohexane']:.5f}", flush=True)

    d = col.design
    dist = fs.streams["DIST"]
    print(f"\nFinal operating point (62 stages, D={d['D']:.1f} mol/s):")
    print(f"  organic reflux (INTERNAL): {d['organic_reflux']:.1f} mol/s, "
          f"x_cyclohexane={d['x_organic']['cyclohexane']:.3f}")
    print(f"  reflux ratio R (output)  : {d['R']:.2f}")
    print(f"  aqueous distillate       : {dist.molar_flow:.1f} mol/s, "
          f"x_water={dist.z['water']:.3f}, x_ethanol={dist.z['ethanol']:.3f}")
    print(f"  bottoms (near-anhydrous) : {fs.streams['BOT'].molar_flow:.1f} mol/s, "
          f"x_ethanol={fs.streams['BOT'].z['ethanol']:.4f}, "
          f"x_water={fs.streams['BOT'].z['water']:.5f}")
    print(f"  reboiler / condenser duty: {d['Q_reboiler'] / 1e3:.1f} / "
          f"{d['Q_condenser'] / 1e3:.1f} kW")


if __name__ == "__main__":
    main()
