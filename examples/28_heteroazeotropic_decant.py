"""Example 28 — heterogeneous-azeotrope (3-phase) decanting.

A heterogeneous azeotrope distillation crosses the distillation boundary in the
*decanter*: the column overhead is the heteroazeotrope (a single vapour), it is
condensed, and it settles into two liquid layers — an organic (entrainer-rich)
layer that is refluxed back to keep the entrainer in the column, and an aqueous
layer drawn as product. That liquid-liquid split is the property package's VLLE
flash; the :class:`~caldyr.unitops.decanter.Decanter` settles it and performs
the reflux-drum split.

This script decants a water / n-butanol heteroazeotrope overhead. It also shows
the new NRTL isoactivity liquid-liquid flash that gives the activity package a
three-phase capability (the cubic EOS gets the structure but overpredicts the
cross-solubility for water/organic systems; an activity model is the right tool
when fit to LLE data).

Run: python examples/28_heteroazeotropic_decant.py
"""
from caldyr.core import Component, Flowsheet
from caldyr.thermo import make_package
from caldyr.unitops import Decanter

P_ATM = 101325.0


def main() -> None:
    # --- the decanter (cubic EOS: robust water/organic structure) ------------
    fs = Flowsheet(components=[Component("water"), Component("n-butanol")],
                   property_package="thermo:PR")
    fs.add(Decanter("DEC", {"T": 320.0, "P": P_ATM,
                            "reflux_layer": "light", "reflux_fraction": 0.7}))
    # a condensed heteroazeotrope overhead (water-rich vapour, now a liquid)
    fs.feed("OVHD", "DEC:in1", T=360.0, P=P_ATM, molar_flow=100.0,
            z={"water": 0.6, "n-butanol": 0.4})
    for port in ("liquid_light", "liquid_heavy", "reflux", "product", "vapor"):
        fs.connect(port.upper(), f"DEC:{port}", None)
    fs.connect("Q", "DEC:duty", None)
    fs.solve()

    dec = fs.units["DEC"]
    refl, prod, aq = (fs.streams["REFLUX"], fs.streams["PRODUCT"],
                      fs.streams["LIQUID_HEAVY"])
    print("Heteroazeotrope decanter (thermo:PR), decant at 320 K:")
    print(f"  organic layer  beta={dec.result['beta_light']:.3f}  "
          f"rho={dec.result['rho_light']:.0f} kg/m^3")
    print(f"  aqueous layer  beta={dec.result['beta_heavy']:.3f}  "
          f"rho={dec.result['rho_heavy']:.0f} kg/m^3")
    print(f"  reflux (organic, 70%): {refl.molar_flow:.1f} mol/s, "
          f"x_butanol={refl.z['n-butanol']:.3f}")
    print(f"  organic product (30%): {prod.molar_flow:.1f} mol/s")
    print(f"  aqueous product:       {aq.molar_flow:.1f} mol/s, "
          f"x_water={aq.z['water']:.3f}")

    # --- the NRTL isoactivity liquid-liquid flash ----------------------------
    # ChemSep NRTL is VLE-fit (no miscibility gap); use illustrative LLE-bearing
    # parameters to show the activity-model VLLE split (organic ~0.46 water /
    # aqueous ~0.99 water — near the experimental 0.51/0.98).
    from thermo import NRTL
    ppn = make_package("thermo:NRTL", ["water", "n-butanol"])
    ppn._flasher.liquid.GibbsExcessModel = NRTL(
        T=298.15, xs=[0.5, 0.5], tau_bs=[[0.0, 1300.0], [1100.0, 0.0]],
        alpha_cs=[[0.0, 0.4], [0.4, 0.0]])
    r = ppn.flash_pt_3p(320.0, P_ATM, {"water": 0.5, "n-butanol": 0.5})
    print("\nNRTL isoactivity VLLE flash (illustrative LLE parameters):")
    print(f"  organic  x_water={r.x_light['water']:.3f}")
    print(f"  aqueous  x_water={r.x_heavy['water']:.3f}")


if __name__ == "__main__":
    main()
