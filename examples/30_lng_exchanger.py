"""Example 30 — multi-stream (LNG) heat exchanger (Hameed 2025 sec. 9.5.2).

The LNG operation solves the heat and material balance for a single plate-fin
core carrying several passes at once — the heart of an LNG liquefaction train,
a cryogenic cold box, or the turbo-expander LPG plant of the reference. Each
pass is listed with a spec (an outlet temperature or a duty, plus a pressure
drop); whether a pass turns out hot or cold is decided by the solution. The
**Weighted** (zone) method builds the hot and cold composite curves with
rigorous PH flashes — capturing any phase change — and reports the minimum
internal approach (MITA) and the required conductance UA.

Here a warm natural-gas feed is cooled by two returning cold streams in one
core. With two outlet temperatures fixed the third pass closes the energy
balance directly; the composite curves then give the MITA and UA. A second run
adds the book's *minimum-approach* specification, freeing a second pass and
solving iteratively until the curves are held exactly 5 K apart.

Run: python examples/30_lng_exchanger.py
"""
from caldyr.core import Component, Flowsheet
from caldyr.economics.costing import cost_equipment
from caldyr.economics.sizing import size_flowsheet
from caldyr.thermo import make_package
from caldyr.unitops import MultiStreamExchanger

P = 4.0e6
COMPS = ["methane", "ethane", "propane", "nitrogen"]
Z = {"methane": 0.88, "ethane": 0.07, "propane": 0.03, "nitrogen": 0.02}


def lng_fs(passes, **extra) -> Flowsheet:
    fs = Flowsheet(components=[Component(c) for c in COMPS],
                   property_package="thermo:PR")
    fs.add(MultiStreamExchanger("LNG", {"passes": passes, **extra}))
    # pass 1: warm feed gas (hot); passes 2 & 3: cold returning gas
    fs.feed("FEED", "LNG:pass1_in", T=310.0, P=P, molar_flow=1.0, z=Z)
    fs.feed("COLD1", "LNG:pass2_in", T=240.0, P=P, molar_flow=0.55, z=Z)
    fs.feed("COLD2", "LNG:pass3_in", T=245.0, P=P, molar_flow=0.55, z=Z)
    for p in ("pass1_out", "pass2_out", "pass3_out"):
        fs.connect(p.upper(), f"LNG:{p}", None)
    return fs


def report(title, passes, **extra):
    fs = lng_fs(passes, **extra)
    rep = fs.solve()
    d = fs.units["LNG"].design
    pp = make_package(fs.property_package, fs.component_ids)
    sizes = size_flowsheet(fs, rep, pp)
    cbm = sum(cost_equipment(s).bare_module for s in sizes)
    print(f"\n{title}")
    print(f"  pass duties (W): "
          f"{', '.join(f'{q:+.0f}' for q in d['pass_duties'])}  "
          f"(Sigma = {sum(d['pass_duties']):+.2e})")
    print(f"  pass outlet T (K): "
          f"{', '.join(f'{t:.1f}' for t in d['pass_T_out'])}")
    print(f"  hot duty {d['hot_duty'] / 1e3:.2f} kW = cold duty "
          f"{d['cold_duty'] / 1e3:.2f} kW")
    print(f"  MITA {d['min_approach']:.2f} K   UA {d['UA'] / 1e3:.2f} kW/K   "
          f"area~{sizes[0].attribute:.2f} m^2   bare-module ${cbm / 1e3:.0f}k")


def main() -> None:
    print("Multi-stream (LNG) exchanger — one plate-fin core, three passes.")
    print("Warm feed gas cooled by two cold returning streams (4 MPa).")

    # (1) Two outlet temperatures fixed -> the third pass closes the balance.
    report("(1) single unknown: pass 3 free, closed by the energy balance",
           [{"T_out": 255.0, "dP": 20000.0}, {"T_out": 300.0}, {}])

    # (2) Add the book's minimum-approach spec -> two passes free, iterative.
    report("(2) min_approach = 5 K spec: passes 2 & 3 free, solved iteratively",
           [{"T_out": 255.0, "dP": 20000.0}, {}, {}], min_approach=5.0)

    # The min_approach run re-balances the two cold passes so the hot and cold
    # composite curves are held exactly 5 K apart at their closest point.


if __name__ == "__main__":
    main()
