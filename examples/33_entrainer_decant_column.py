"""Example 33 — anhydrous ethanol via an integrated DECANTING-CONDENSER column.

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
in full, and the aqueous layer is the distillate (spec = ``distillate_rate``). The
reflux ratio and both duties are recovered outputs.

Removing **all** the feed water (anhydrous bottoms) still needs a large aqueous
draw — water leaves only in the ethanol-rich aqueous layer — so the operating
point sits deep in the large-distillate regime. We reach it by **warm-started
distillate-rate continuation**: a cold solve at a small, easy D, then ramp D up
(each step re-converges in a handful of Newton iterations from the previous
profile). At D≈18 the bottoms is ~98.7 % ethanol / ~0.8 % water / ~0.5 %
cyclohexane — the entrainer is recirculating internally and the column is in the
anhydrous regime the external decanter could not converge.

Run: python examples/33_entrainer_decant_column.py
"""
from caldyr.core import Component, Flowsheet
from caldyr.unitops import RigorousColumn

P_ATM = 101325.0


def main() -> None:
    col = RigorousColumn("T100", {
        "n_stages": 30,
        "feeds": [{"stage": 3}, {"stage": 14}],     # cyclohexane top, feed central
        "reflux_ratio": 3.0,                        # seed organic reflux ratio
        "distillate_rate": 4.2,                     # ramped up below
        "method": "naphtali_sandholm",
        "reboiled": True,
        "decant_condenser": True,
        "condenser_T": 305.0,                       # decant the overhead at 32 C
        "reflux_layer": "organic",                  # reflux the cyclohexane layer
        "max_iter": 120,
    })

    # Warm-started distillate-rate continuation toward the anhydrous regime.
    print("Integrated decanting-condenser column (ethanol/water/cyclohexane):")
    print(f"  {'D':>5} {'makeup':>7} {'iters':>6} {'Eres':>9}   "
          f"{'bottoms ethanol / water / cyclohexane':<40}")
    for D, solv in [(4.2, 8.0), (6.0, 7.0), (8.0, 6.0), (10.0, 5.0),
                    (12.0, 4.0), (14.0, 3.5), (16.0, 3.2), (18.0, 3.0)]:
        col.params["distillate_rate"] = D
        fs = Flowsheet(
            components=[Component("ethanol"), Component("water"),
                        Component("cyclohexane")],
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
        d = col.design
        bot = fs.streams["BOT"]
        print(f"  {D:5.1f} {solv:7.1f} {d['iterations']:6d} "
              f"{d['energy_residual_rel']:9.1e}   "
              f"{bot.z['ethanol']:.4f} / {bot.z['water']:.4f} / "
              f"{bot.z['cyclohexane']:.4f}")

    d = col.design
    dist = fs.streams["DIST"]
    print(f"\nFinal operating point (D={d['D']:.1f} mol/s):")
    print(f"  organic reflux (INTERNAL): {d['organic_reflux']:.1f} mol/s, "
          f"x_cyclohexane={d['x_organic']['cyclohexane']:.3f}")
    print(f"  reflux ratio R (output)  : {d['R']:.2f}")
    print(f"  aqueous distillate       : {dist.molar_flow:.1f} mol/s, "
          f"x_water={dist.z['water']:.3f}, x_ethanol={dist.z['ethanol']:.3f}")
    print(f"  bottoms (anhydrous EtOH) : {fs.streams['BOT'].molar_flow:.1f} mol/s, "
          f"x_ethanol={fs.streams['BOT'].z['ethanol']:.4f}")
    print(f"  reboiler / condenser duty: {d['Q_reboiler'] / 1e3:.1f} / "
          f"{d['Q_condenser'] / 1e3:.1f} kW")


if __name__ == "__main__":
    main()
