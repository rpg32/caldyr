"""M18 demo: reactive distillation (Hameed §9.5.3).

A reactive distillation column runs a reaction on a band of trays while it
separates the products at the same time. Continuously pulling a product out of
the reactive zone shifts a reversible equilibrium toward more conversion (Le
Chatelier) — higher yield in one vessel, the attraction of the technique.

``RigorousColumn`` takes kinetic reactions (forward and, for an equilibrium-
limited system, reversible — see ``reactions`` / ``tray_holdup``) on named
stages: the per-stage generation enters the MESH component balances as an extra
source and the heat of reaction is carried by the formation-inclusive stage
enthalpies.

Part 1 — **methyl-acetate synthesis** (the book's example): methanol + acetic
acid ⇌ methyl acetate + water, on stages 5-10 of a 15-stage column. The book uses
HYSYS/Wilson; caldyr's activity-model analogue is ``thermo:NRTL``. The reaction
makes the light ester, which concentrates in the distillate.

Part 2 — **toluene disproportionation** (2 toluene ⇌ benzene + xylene), a clean
all-aromatic cross-check on ``thermo:PR`` that reaches ~60 % conversion with
benzene overhead and xylene as bottoms.

    python examples/27_reactive_distillation.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.unitops import RigorousColumn  # noqa: E402


def methyl_acetate(holdup: float) -> Flowsheet:
    rxn = {"stoich": {"methanol": -1, "acetic acid": -1,
                      "methyl acetate": 1, "water": 1},
           "key": "methanol", "k0": 3e-8, "Ea": 0.0,
           "orders": {"methanol": 1, "acetic acid": 1},
           "k0_rev": 6e-9, "Ea_rev": 0.0,
           "orders_rev": {"methyl acetate": 1, "water": 1},
           "stages": [5, 10]}
    params = dict(n_stages=15, feed_stage=10, reflux_ratio=5.0,
                  distillate_rate=5.556, P=90000.0, dP_stage=500.0)
    if holdup > 0:
        params["reactions"] = [rxn]
        params["tray_holdup"] = holdup
    comps = ["methanol", "acetic acid", "methyl acetate", "water"]
    fs = Flowsheet(components=[Component(c) for c in comps],
                   property_package="thermo:NRTL")
    fs.add(RigorousColumn("COL", params))
    fs.feed("F", "COL:in1", T=348.15, P=101325.0, molar_flow=12.5,
            z={"methanol": 0.4, "acetic acid": 0.4,
               "methyl acetate": 0.1, "water": 0.1})
    for n, port in [("D", "distillate"), ("B", "bottoms"),
                    ("qc", "condenser_duty"), ("qr", "reboiler_duty")]:
        fs.connect(n, "COL:" + port, None)
    return fs


def toluene(holdup: float) -> Flowsheet:
    rxn = {"stoich": {"toluene": -2, "benzene": 1, "p-xylene": 1},
           "key": "toluene", "k0": 3e-8, "Ea": 0.0, "orders": {"toluene": 2},
           "k0_rev": 7.5e-8, "Ea_rev": 0.0,
           "orders_rev": {"benzene": 1, "p-xylene": 1}, "stages": [4, 11]}
    params = dict(n_stages=15, feed_stage=8, reflux_ratio=3.0,
                  distillate_to_feed=0.5)
    if holdup > 0:
        params["reactions"] = [rxn]
        params["tray_holdup"] = holdup
    comps = ["benzene", "toluene", "p-xylene"]
    fs = Flowsheet(components=[Component(c) for c in comps],
                   property_package="thermo:PR")
    fs.add(RigorousColumn("COL", params))
    fs.feed("F", "COL:in1", T=400.0, P=120000.0, molar_flow=10.0,
            z={"benzene": 0.05, "toluene": 0.90, "p-xylene": 0.05})
    for n, port in [("D", "distillate"), ("B", "bottoms"),
                    ("qc", "condenser_duty"), ("qr", "reboiler_duty")]:
        fs.connect(n, "COL:" + port, None)
    return fs


def main() -> None:
    print("=== Part 1: methyl-acetate synthesis (book §9.5.3, NRTL) ===")
    for holdup in (0.0, 0.3):
        fs = methyl_acetate(holdup)
        fs.solve()
        d, b = fs.streams["D"], fs.streams["B"]
        meoh_out = d.molar_flow * d.z["methanol"] + b.molar_flow * b.z["methanol"]
        conv = (12.5 * 0.4 - meoh_out) / (12.5 * 0.4)
        tag = "no reaction" if holdup == 0 else f"holdup {holdup} m3"
        print(f"  {tag:14s}: MeOH conversion {conv * 100:5.1f} %, "
              f"distillate x_MeAc {d.z['methyl acetate']:.3f}, "
              f"bottoms x_MeAc {b.z['methyl acetate']:.3f}")

    print("\n=== Part 2: toluene disproportionation (PR) ===")
    for holdup in (0.0, 0.5):
        fs = toluene(holdup)
        fs.solve()
        d, b = fs.streams["D"], fs.streams["B"]
        tol_out = d.molar_flow * d.z["toluene"] + b.molar_flow * b.z["toluene"]
        conv = (10.0 * 0.9 - tol_out) / (10.0 * 0.9)
        tag = "no reaction" if holdup == 0 else f"holdup {holdup} m3"
        print(f"  {tag:14s}: toluene conversion {conv * 100:5.1f} %, "
              f"distillate x_benzene {d.z['benzene']:.3f}, "
              f"bottoms x_xylene {b.z['p-xylene']:.3f}")


if __name__ == "__main__":
    main()
